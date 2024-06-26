#' Convert a network dynamic object into a cumulative edgelist
netdyn2el_cuml <- function(net) {
  as.data.frame(net) |>
    dplyr::select(start = onset, stop = terminus, head, tail) |>
    dplyr::mutate(stop = stop - 1)
}

#' @title Calculate the Forward Reachable Path over a Time Series
#'
#' @description This function calculates the Forward Reachable Path (FRP) of all
#'              the nodes in a network over a time series. It is much faster
#'              than iterating \code{tsna::tPath} over all nodes.
#'
#' @param el_cuml a cumulative edgelist object. That is a data.frame with at
#'   least columns: head, tail, start and stop. Start and stop are inclusive.
#' @param from_step the beginning of the time period.
#' @param to_step the end of the time period.
#'
#' @return
#' A matrix of list with \code{n_nodes} rows and \code{n_steps + 1}
#' columns, containing the nodes adding to the FRP of each node at each
#' time step. The first columns contains only the nodes themselves, and each
#' subsequent columns contains the list of the nodes that are added to the FRP
#' at time step \code{from + colnum + 1}
#'
#' @details
#' See the examples for how to recover the full FRP of a node between
#' \code{from} and any time step up to \code{end}.
#'
#' @section Time and Memory Use:
#' This function may be used to efficiently calculate all FRPs over many time
#' steps. For more limited calculations, see \code{tsna::tPath}. This function
#' takes 3 to 20 minutes on a network of 1e4 nodes over 260 time steps.
#'
#' @section Displaying Progress:
#' This function is using the
#' \href{https://progressr.futureverse.org/articles/progressr-intro.html}{progressr package}
#' to display its progression. Use
#' \code{progressr::with_progress({frp_parts <- get_all_frp(net, from = 1, to = 260)})}
#' to display the progress bar. Or see the
#' \href{https://progressr.futureverse.org/articles/progressr-intro.html}{progressr package}
#' for more information and customization.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate all the FRPs from step 100 to 260
#' from_ts <- 100
#' to_ts <- 260
#'
#' frp_parts <- get_all_frp(el_cuml, from_step = from_ts, to_stop = to_ts)
#'
#' # Get the FRP of node 10 from step 100 (from_ts) to 150
#' # note how from_ts does not appear as we can only calculate from it once
#' # the function has been run.
#' frp_10 <- unlist(frp_parts[10, seq_len(50 + 1)])
#'
#' # Get the length of the FRPs for each node at each timestep
#' frp_parts_length <- apply(frp_parts, c(1, 2), function(x) length(x[[1]]))
#' frp_lengths <- t(apply(frp_parts_length , 1,cumsum))
#'
#' # testing the results against tPath
#' n_max <- 500
#' n <- 0
#' while(n < n_max) {
#'   v_int <- sample(n_nodes, 1)
#'   ts <- sample(n_steps, 1)
#'
#'   # get the FRP using tPath
#'   tp <- tsna::tPath(net, v = v_int,
#'                     start = from_ts, end = from_ts + ts,
#'                     direction = "fwd")
#'   frp_tp <- which(tp$tdist < Inf)
#'
#'   # get the FRP using this function
#'   frp_my <- unlist(frp_parts[v_int, seq_len(ts + 1)])
#'
#'   if (!setequal(frp_tp, frp_my))
#'     stop("missmatch in node: ", v_int, "; for ts = ", ts)
#'   n <- n + 1
#'   print(n)
#' }
#'
#' }
get_all_frp <- function(el_cuml, from_step, to_step) {
  n_steps <- to_step - from_step + 1
  n_nodes <- max(c(el_cuml$head, el_cuml$tail))

  # nolint start
  df_net <- el_cuml |>
    dplyr::mutate(stop = ifelse(is.na(stop), Inf, stop)) |> # current edges never ends
    dplyr::filter(start <= to_step, stop >= from_step) |> # remove edges before and after analysis period
    dplyr::select(start, stop, head, tail)

  # nolint end
  # the initial FRP contains only the vertex itself
  frp_cur <- as.list(seq_len(n_nodes))
  frp_parts <- matrix(list(numeric(0)), ncol = n_steps + 1, nrow = n_nodes)
  frp_parts[, 1] <- frp_cur

  p <- progressr::progressor(n_steps)
  for (t in seq_len(n_steps)) {
    p()
    cur_step <- t + from_step - 1

    # creation of a `connection` list of vectors
    # for each vertex 1:n_nodes we get a vector of the vertices it connects to
    connected <- vector(mode = "list", length = n_nodes)
    # nolint start
    # IN EL_CUML: duration is [start, stop] (inclusive)
    # all current edges are needed (not only start). This is because we use the
    # connection to each node to calc the FRP to the other ones. Only start
    # would be enough for single FRP though
    el_t <- dplyr::filter(df_net, start <= cur_step, stop >= cur_step)
    # nolint end
    for (i in seq_len(nrow(el_t))) {
      e_head <- el_t$head[i]
      e_tail <- el_t$tail[i]
      connected[[e_head]] <- c(connected[[e_head]], e_tail)
      connected[[e_tail]] <- c(connected[[e_tail]], e_head)
    }

    # PERF: bottleneck is here
    # frp_v is the current frp for vertex v at timestep t - 1
    # we add to it all the nodes that have edges at timestep t with any of the
    # nodes in the FRP
    # the while loop is to include the nodes that are connected to the FRP through
    # a node added this step

    # use hashmaps? rust?
    frp_new <- lapply(
      frp_cur,
      function(frp_v) {
        only_new <- numeric(0)
        new <- frp_v
        while (length(new) > 0 & length(frp_v) < n_nodes) {
          new <- unlist(connected[new])
          new <- setdiff(new, frp_v)
          frp_v <- c(frp_v, new)
          only_new <- c(only_new, new)
        }
        only_new
      }
    )

    # could speed up here with "true arrays" (growable, pre-alloc)
    frp_cur <- Map(c, frp_cur, frp_new)
    frp_parts[, t + 1] <- frp_new
  }
  return(frp_parts)
}

# re-create the frp_cur per node at each step
# slower but less memory intensive
get_all_frp_nocur <- function(el_cuml, from_step, to_step) {
  n_steps <- to_step - from_step + 1
  n_nodes <- max(c(el_cuml$head, el_cuml$tail))

  # nolint start
  df_net <- el_cuml |>
    dplyr::mutate(stop = ifelse(is.na(stop), Inf, stop)) |> # current edges never ends
    dplyr::filter(start <= to_step, stop >= from_step) |> # remove edges before and after analysis period
    dplyr::select(start, stop, head, tail)

  # nolint end
  # the initial FRP contains only the vertex itself
  frp_parts <- matrix(list(numeric(0)), ncol = n_steps + 1, nrow = n_nodes)
  frp_parts[, 1] <- as.list(seq_len(n_nodes))

  p <- progressr::progressor(n_steps)
  for (t in seq_len(n_steps)) {
    p()
    cur_step <- t + from_step - 1

    # creation of a `connection` list of vectors
    # for each vertex 1:n_nodes we get a vector of the vertices it connects to
    connected <- vector(mode = "list", length = n_nodes)
    # nolint start
    # IN EL_CUML: duration is [start, stop] (inclusive)
    # all current edges are needed (not only start). This is because we use the
    # connection to each node to calc the FRP to the other ones. Only start
    # would be enough for single FRP though
    el_t <- dplyr::filter(df_net, start <= cur_step, stop >= cur_step)
    # nolint end
    for (i in seq_len(nrow(el_t))) {
      e_head <- el_t$head[i]
      e_tail <- el_t$tail[i]
      connected[[e_head]] <- c(connected[[e_head]], e_tail)
      connected[[e_tail]] <- c(connected[[e_tail]], e_head)
    }

    # PERF: bottleneck is here
    # frp_v is the current frp for vertex v at timestep t - 1
    # we add to it all the nodes that have edges at timestep t with any of the
    # nodes in the FRP
    # the while loop is to include the nodes that are connected to the FRP through
    # a node added this step

    # use hashmaps? rust?
    frp_new <- lapply(
      seq_len(n_nodes),
      function(i) {
        frp_v <- unlist(frp_parts[i, 1:t])
        only_new <- numeric(0)
        new <- frp_v
        while (length(new) > 0 & length(frp_v) < n_nodes) {
          new <- unlist(connected[new])
          new <- setdiff(new, frp_v)
          frp_v <- c(frp_v, new)
          only_new <- c(only_new, new)
        }
        only_new
      }
    )

    # could speed up here with "true arrays" (growable, pre-alloc)
    frp_parts[, t + 1] <- frp_new
  }
  return(frp_parts)
}

get_all_frp_old <- function(net, from = 1, to = Inf) {
  last_obs <- length(net$gal$net.obs.period$observations)
  to <- if (to > last_obs) last_obs else to
  n_steps <- to - from + 1
  n_nodes <- net$gal$n

  # nolint start
  onset <- terminus <- head <- tail <- NULL
  df_net <- dplyr::select(
    as.data.frame(net),
    onset, terminus, head, tail
  )
  # nolint end
  # the initial FRP contains only the vertex itself
  frp_cur <- as.list(seq_len(n_nodes))
  frp_parts <- matrix(list(numeric(0)), ncol = n_steps + 1, nrow = n_nodes)
  frp_parts[, 1] <- frp_cur

  p <- progressr::progressor(n_steps)
  for (t in seq_len(n_steps)) {
    p()
    cur_step <- t + from - 1

    # creation of a `connection` list of vectors
    # for each vertex 1:n_nodes we get a vector of the vertices it connects to
    connected <- vector(mode = "list", length = n_nodes)
    # nolint start
    el_t <- dplyr::filter(df_net, onset <= cur_step, terminus > cur_step)
    el_t <- dplyr::select(el_t, head, tail)
    # nolint end
    for (i in seq_len(nrow(el_t))) {
      e_head <- el_t$head[i]
      e_tail <- el_t$tail[i]
      connected[[e_head]] <- c(connected[[e_head]], e_tail)
      connected[[e_tail]] <- c(connected[[e_tail]], e_head)
    }

    # PERF: bottleneck is here
    # frp_v is the current frp for vertex v at timestep t - 1
    # we add to it all the nodes that have edges at timestep t with any of the
    # nodes in the FRP
    # the while loop is to include the nodes that are connected to the FRP through
    # a node added this step
    frp_new <- lapply(
      frp_cur,
      function(frp_v) {
        only_new <- numeric(0)
        new <- frp_v
        while (length(new) > 0 & length(frp_v) < n_nodes) {
          new <- unlist(connected[new])
          new <- setdiff(new, frp_v)
          frp_v <- c(frp_v, new)
          only_new <- c(only_new, new)
        }
        only_new
      }
    )

    frp_cur <- Map(c, frp_cur, frp_new)
    frp_parts[, t + 1] <- frp_new
  }

  return(frp_parts)
}
