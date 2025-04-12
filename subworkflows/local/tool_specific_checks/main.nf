workflow TOOL_SPECIFIC_CHECKS {
    main:
    def checkPathParamList = []

    // Conditionally add tools specific parameters (only examples for now)
    if (params.tools?.split(',')?.contains('snpeff') || params.tools?.split(',')?.contains('merge')) {
        checkPathParamList.add(params.snpeff_cache)
    }
    if (params.tools?.split(',')?.contains('vep') || params.tools?.split(',')?.contains('merge')) {
        checkPathParamList.add(params.vep_cache)
    }

    if (checkPathParamList) {
        paths_to_check_channel = Channel.fromList(checkPathParamList)

        missing_paths = paths_to_check_channel.map { path_item ->
            def file_path = path_item instanceof Path ? path_item : file(path_item)
            return [file_path, file_path.exists()]
        }
        .filter { !it[1] }
        .collect()
        .subscribe { missing ->
            if (missing.size() > 0) {
               def error_string = "The following required paths do not exist:\n" +
                                   missing.collect { "  - ${it[0]}" }.join('\n')
               error(error_string)
            }
        }
    }
}
