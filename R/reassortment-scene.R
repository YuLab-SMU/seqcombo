prepare_reassortment_scene <- function(virus_info, flow_info = NULL, asp = 1,
                                       v_shape = "ellipse") {
    validate_virus_info(virus_info, require_coordinates = TRUE)
    validate_segment_color(virus_info)
    validate_segment_name(virus_info)

    if (!is.null(flow_info)) {
        validate_flow_info(flow_info, virus_info$id)
    }

    ASP <- asp_(virus_info, asp)
    virus_info <- set_virus_size(virus_info, ASP, v_shape)
    capsule_data <- get_capsule_data(virus_info, ASP, v_shape)

    label_data <- NULL
    if (all(c("label", "label_position") %in% colnames(virus_info))) {
        label_data <- generate_label_data(virus_info, capsule_data)
    }

    flow_data <- NULL
    if (!is.null(flow_info)) {
        flow_data <- route_reassortment_edges(virus_info, flow_info, capsule_data, ASP)
        extra_flow_col <- setdiff(colnames(flow_info), c("from", "to"))
        if (length(extra_flow_col) > 0) {
            flow_data <- cbind(flow_data, flow_info[, extra_flow_col, drop = FALSE])
        }
    }

    list(
        virus_info = virus_info,
        flow_info = flow_info,
        asp = ASP,
        v_shape = v_shape,
        capsule_data = capsule_data,
        label_data = label_data,
        flow_data = flow_data
    )
}
