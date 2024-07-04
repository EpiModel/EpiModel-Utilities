forward.reachable <- function(nd, v, start = NULL, end = NULL,
                              per.step.depth = Inf) {
  if (!is.networkDynamic(nd)) {
    stop("the first argument to forward.reachble must be a networkDynamic object")
  }
  if (missing(v) || !is.numeric(v)) {
    stop("v argument to forward.reachable must be a vector of valid numeric vertex ids")
  }
  if (max(v) > network.size(nd) | min(v) < 1) {
    stop("v argument to forward.reachable must be a vector of numeric vertex ids within the range of the network size")
  }

  # set the interval to be whatever the observed changes are
  times <- get.change.times(
    nd,
    vertex.attribute.activity = FALSE,
    edge.attribute.activity = FALSE,
    network.attribute.activity = FALSE
  )

  if (length(times) == 0) {
    times <- c(0, Inf)
  }
  if (is.null(start)) {
    # start<-min(times)
    start <- -Inf
  }
  if (is.null(end)) {
    # end<-max(times)
    end <- Inf
  }
  # trim times to desired range, making sure to include start and end
  times <- unique(c(start, times[times >= start]))
  times <- unique(c(times[times <= end], end))

  distance <- rep(Inf, network.size(nd))
  distance[v] <- times[1]

  # TODO: could probably skip all times earlier that the active times in v?
  reached <- v
  for (t in 1:(length(times) - 1)) {
    # BFS to depth rate
    new <- reached
    # how long until next change?
    duration <- times[t + 1] - times[t]

    # remove any in the set we've already visited
    if (duration > 0) {
      d <- 1 # we are assuming all geodesic steps count as 1, harder if we calc per edge..
      # keep searching until we reach bounds or run out of verts to find
      # also stop if we find all the vertices
      while (d <= per.step.depth * duration & length(reached) < network.size(nd)) {
        ngs <- unlist(unique(sapply(new, function(i) {
          get.neighborhood.active(nd, v = i, at = times[t], type = "out")
        })))
        new <- setdiff(ngs, reached)
        if (length(new) == 0) {
          break # no more verts to find
        }
        distance[new] <- times[t]
        reached <- c(reached, new)
        d <- d + 1
      }
    }
  }
  return(reached)
}

get_frp <- function(el_cuml, node, from_step, to_step) {
  # nolint start
  df_net <- el_cuml |>
    dplyr::mutate(stop = ifelse(is.na(stop), Inf, stop)) |> # current edges never ends
    dplyr::filter(start <= to_step, stop >= from_step) |> # remove edges before and after analysis period
    dplyr::mutate(start = ifelse(start < from_step, from_step, start)) |> # set older edges to start at beginning of analysis
    dplyr::select(start, stop, head, tail)

  # nolint end
  change_times <- sort(unique(df_net$start))

  # the initial FRP contains only the vertex itself
  frp_cur <- as.list(seq_len(n_nodes))
  frp_parts <- matrix(list(numeric(0)), ncol = n_steps + 1, nrow = n_nodes)
  frp_parts[, 1] <- frp_cur

  p <- progressr::progressor(length(change_times))
  for (cur_step in change_times) {
    p()
    t <- cur_step - from_step + 1 # the slot in the frp_parts

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

    only_new <- numeric(0)
    new <- frp_v
    while (length(new) > 0 && length(frp_v) < n_nodes) {
      new <- unlist(connected[new])
      new <- setdiff(new, frp_v)
      frp_v <- c(frp_v, new)
      only_new <- c(only_new, new)
    }
    only_new

    # could speed up here with "true arrays" (growable, pre-alloc)
    frp_cur <- Map(c, frp_cur, frp_new)
    frp_parts[, t + 1] <- frp_new
  }
  return(frp_parts)
}
