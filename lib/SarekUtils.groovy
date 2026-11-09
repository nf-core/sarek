class SarekUtils {
    // Case-insensitive membership test for comma-separated --tools / --skip_tools values
    static boolean hasTool(String list, String tool) {
        list?.toLowerCase()?.tokenize(',')?.contains(tool.toLowerCase()) ?: false
    }
}
