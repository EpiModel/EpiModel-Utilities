#' Convert a network dynamic object into a cumulative edgelist
netdyn2el_cuml <- function(net) {
  as.data.frame(net) |>
    dplyr::select(start = onset, stop = terminus, head, tail) |>
    dplyr::mutate(stop = stop - 1)
}

# dedplicate a cumulative edgelist by combining overllapping edges into one
# with the full duration of the overllapping ones
dedup_el_cuml <- function(el_all) {
  ea <- el_all |>
    dplyr::group_by(head, tail) |>
    dplyr::mutate(n = n())

  e_unique <- ea |>
    dplyr::filter(n == 1) |>
    dplyr::select(-n)

  e_dup <- ea |>
    dplyr::filter(n > 1) |>
    dplyr::select(-n)

  e_dup <- e_dup |>
    dplyr::arrange(head, tail, start, stop) |>
    dplyr::group_by(head, tail)

  e_dedup <- e_dup |>
    dplyr::mutate(
      lstart= dplyr::lag(start),
      lstop = dplyr::lag(stop),
      overlap = !is.na(lstop) & !is.na(lstart) & start <= lstop,
      stop = ifelse(overlap, max(stop, lstop, na.rm = TRUE), stop),
      start = ifelse(overlap, min(start, lstart, na.rm = TRUE), start)
      ) |>
    dplyr::select(-c(lstart, lstop, overlap)) |>
    dplyr::ungroup() |>
    unique()

  dplyr::bind_rows(e_unique, e_dedup)
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
#' @param nodes the subset of nodes to calculate the FRP for. (default = NULL,
#'        all nodes)
#'
#' @return
#' A list of FRP for each of the nodes of interest
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
#' @section Number of Nodes:
#' This codes does not know the total number of node on the network and assumes
#' that the highest ID recorded correspond to the last node.
#' We can therefore arrive to a situation where there is elements in the output
#' than node in the network if the last N nodes (by ID) are never connected.
#' And therefore are not recorded in the cumulative edgelist.
#' So the FRP for the nodes not present in the output is always 1 (themselves).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate all the FRPs from step 100 to 260
#' from_ts <- 100
#' to_ts <- 260
#'
#' frps <- get_all_frp(el_cuml, from_step = from_ts, to_stop = to_ts)
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
#'   frp_my <- frps[[v_int]]
#'
#'   if (!setequal(frp_tp, frp_my))
#'     stop("missmatch in node: ", v_int, "; for ts = ", ts)
#'   n <- n + 1
#'   print(n)
#' }
#'
#' }
get_all_frp <- function(el_cuml, from_step, to_step, nodes = NULL) {
  n_nodes <- max(c(el_cuml$head, el_cuml$tail))

  # nolint start
  df_net <- el_cuml |>
    dplyr::mutate(stop = ifelse(is.na(stop), Inf, stop)) |> # current edges never ends
    dplyr::filter(start <= to_step, stop >= from_step) |> # remove edges before and after analysis period
    dplyr::mutate(start = ifelse(start < from_step, from_step, start)) |> # set older edges to start at beginning of analysis
    dplyr::select(start, stop, head, tail)
  # nolint end
  change_times <- sort(unique(df_net$start))

  if (is.null(nodes))
    nodes <- seq_len(n_nodes)

  # the initial FRP contains only the vertex itself
  frp_cur <- as.list(nodes)
  names(frp_cur) <- paste0("node_", nodes)

  p <- progressr::progressor(length(change_times))
  for (cur_step in change_times) {
    p()
    t <- cur_step - from_step + 1 # the slot in the frp_parts

    # IN EL_CUML: duration is [start, stop] (inclusive)
    # all current edges are needed (not only start). This is because getting
    # connected to a node A means that we are also indirectly connected to its
    # connections, even the ones that started previously.

    # nolint start
    el_t <- dplyr::filter(df_net, start <= cur_step, stop >= cur_step)
    # nolint end

    # creation of a `connection` list of vectors
    # for each vertex 1:n_nodes we get a vector of the vertices it connects to
    connected <- get_connected(el_t, n_nodes)

    # PERF: bottleneck is here
    # frp_v is the current frp for vertex v at timestep t - 1
    # we add to it all the nodes that have edges at timestep t with any of the
    # nodes in the FRP
    # the while loop is to include the nodes that are connected to the FRP
    # through a node added this step

    frp_cur <- lapply(frp_cur, get_subnet,
                      connected = connected, n_nodes = n_nodes)
  }
  return(frp_cur)
}

get_frp_lengths <- function(el_cuml, from_step, to_step, nodes = NULL) {
  frps <- get_all_frp(el_cuml, from_step, to_step, nodes)
  vapply(frps, length, numeric(1))
}

# `nodes` are connecting the parts of the subnet if disjointed
# e.g. `nodes` are the previous FRP
get_subnet <- function(connected, nodes, n_nodes) {
  new <- nodes
  subnet <- nodes
  while (length(new) > 0 && length(subnet) < n_nodes) {
    new <- unlist(connected[new])
    new <- setdiff(new, subnet)
    subnet <- unique(c(subnet, new))
  }
  subnet
}

# Make a `connection` fast acces list from an edgelist
get_connected <- function(el, n_nodes) {
  connected <- vector(mode = "list", length = n_nodes)
  for (i in seq_len(nrow(el))) {
    e_head <- el$head[i]
    e_tail <- el$tail[i]
    connected[[e_head]] <- c(connected[[e_head]], e_tail)
    connected[[e_tail]] <- c(connected[[e_tail]], e_head)
  }
  connected
}

