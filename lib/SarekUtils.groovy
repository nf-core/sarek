import static nextflow.Nextflow.file
import nextflow.Channel

class SarekUtils {

  // Check parameter existence
  static def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
      println("Unknown parameter: ${it}")
      return false
    }
    return true
  }

  // Compare each parameter with a list of parameters
  static def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
  }

  // Loop through all the references files to check their existence
  static def checkReferenceMap(referenceMap) {
    referenceMap.every {
      referenceFile, fileToCheck ->
      SarekUtils.checkRefExistence(referenceFile, fileToCheck)
    }
  }

  // Loop through all the references files to check their existence
  static def checkRefExistence(referenceFile, fileToCheck) {
    if (fileToCheck instanceof List) return fileToCheck.every{ SarekUtils.checkRefExistence(referenceFile, it) }
    def f = file(fileToCheck)
    // this is an expanded wildcard: we can assume all files exist
    if (f instanceof List && f.size() > 0) return true
    else if (!f.exists()) {
			println  "Missing references: ${referenceFile} ${fileToCheck}"
      return false
    }
    return true
  }

}
