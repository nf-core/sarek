
// Global function
// Check the parameters tools or skip_tools, then compare it against the provided tool
// Returns true/false based on whether 'tool' is found in 'parameter'
def checkInParam(parameter, checkValue) {
    if (!parameter){
        false
    } else {
        def tokenized_parameter = parameter.tokenize(',')
        switch (checkValue) {
            // If checkValue is a list check if any appear in tokenized parameter
            case checkValue instanceof List:
                checkValue.any{ it.toLowerCase() in parameter.tokenize(',') }
            // If checkValue is a string check it appears in parameter
            case checkValue instanceof String:
                checkValue.toLowerCase() in parameter.tokenize(',')
            default:
                false
        }
    }
}