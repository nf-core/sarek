//
// Subworkflow with utility functions specific to the nf-core pipeline template
//

/*
========================================================================================
    SUBWORKFLOW DEFINITION
========================================================================================
*/

workflow UTILS_NFCORE_PIPELINE {

    main:
    checkConfigProvided()

    emit:
    success = true

}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
//  Warn if a -profile or Nextflow config has not been provided to run the pipeline
//
def checkConfigProvided() {
    if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
        log.warn "[$workflow.manifest.name] You are attempting to run the pipeline without any custom configuration!\n\n" +
            "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" +
            "   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`\n" +
            "   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`\n" +
            "   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`\n\n" +
            "Please refer to the quick start section and usage docs for the pipeline.\n "
    }
}

//
// Citation string for pipeline
//
def workflowCitation() {
    return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
        "* The pipeline\n" +
        "  ${workflow.manifest.doi}\n\n" +
        "* The nf-core framework\n" +
        "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
        "* Software dependencies\n" +
        "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
}

//
// nf-core logo
//
def nfCoreLogo(workflow_version) {
    Map colors = logColours()
    String.format(
        """\n
        ${dashedLine()}
                                                ${colors.green},--.${colors.black}/${colors.green},-.${colors.reset}
        ${colors.blue}        ___     __   __   __   ___     ${colors.green}/,-._.--~\'${colors.reset}
        ${colors.blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${colors.yellow}}  {${colors.reset}
        ${colors.blue}  | \\| |       \\__, \\__/ |  \\ |___     ${colors.green}\\`-._,-`-,${colors.reset}
                                                ${colors.green}`._,._,\'${colors.reset}
        ${colors.purple}  ${workflow.manifest.name} ${workflow_version}${colors.reset}
        ${dashedLine()}
        """.stripIndent()
    )
}

//
// Return dashed line
//
def dashedLine() {
    Map colors = logColours()
    return "-${colors.dim}----------------------------------------------------${colors.reset}-"
}

//
// ANSII colours used for terminal logging
//
def logColours() {
    Map colorcodes = [:]

    // Reset / Meta
    colorcodes['reset']      = params.monochrome_logs ? '' : "\033[0m"
    colorcodes['bold']       = params.monochrome_logs ? '' : "\033[1m"
    colorcodes['dim']        = params.monochrome_logs ? '' : "\033[2m"
    colorcodes['underlined'] = params.monochrome_logs ? '' : "\033[4m"
    colorcodes['blink']      = params.monochrome_logs ? '' : "\033[5m"
    colorcodes['reverse']    = params.monochrome_logs ? '' : "\033[7m"
    colorcodes['hidden']     = params.monochrome_logs ? '' : "\033[8m"

    // Regular Colors
    colorcodes['black']      = params.monochrome_logs ? '' : "\033[0;30m"
    colorcodes['red']        = params.monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']      = params.monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow']     = params.monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']       = params.monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple']     = params.monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']       = params.monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']      = params.monochrome_logs ? '' : "\033[0;37m"

    // Bold
    colorcodes['bblack']     = params.monochrome_logs ? '' : "\033[1;30m"
    colorcodes['bred']       = params.monochrome_logs ? '' : "\033[1;31m"
    colorcodes['bgreen']     = params.monochrome_logs ? '' : "\033[1;32m"
    colorcodes['byellow']    = params.monochrome_logs ? '' : "\033[1;33m"
    colorcodes['bblue']      = params.monochrome_logs ? '' : "\033[1;34m"
    colorcodes['bpurple']    = params.monochrome_logs ? '' : "\033[1;35m"
    colorcodes['bcyan']      = params.monochrome_logs ? '' : "\033[1;36m"
    colorcodes['bwhite']     = params.monochrome_logs ? '' : "\033[1;37m"

    // Underline
    colorcodes['ublack']     = params.monochrome_logs ? '' : "\033[4;30m"
    colorcodes['ured']       = params.monochrome_logs ? '' : "\033[4;31m"
    colorcodes['ugreen']     = params.monochrome_logs ? '' : "\033[4;32m"
    colorcodes['uyellow']    = params.monochrome_logs ? '' : "\033[4;33m"
    colorcodes['ublue']      = params.monochrome_logs ? '' : "\033[4;34m"
    colorcodes['upurple']    = params.monochrome_logs ? '' : "\033[4;35m"
    colorcodes['ucyan']      = params.monochrome_logs ? '' : "\033[4;36m"
    colorcodes['uwhite']     = params.monochrome_logs ? '' : "\033[4;37m"

    // High Intensity
    colorcodes['iblack']     = params.monochrome_logs ? '' : "\033[0;90m"
    colorcodes['ired']       = params.monochrome_logs ? '' : "\033[0;91m"
    colorcodes['igreen']     = params.monochrome_logs ? '' : "\033[0;92m"
    colorcodes['iyellow']    = params.monochrome_logs ? '' : "\033[0;93m"
    colorcodes['iblue']      = params.monochrome_logs ? '' : "\033[0;94m"
    colorcodes['ipurple']    = params.monochrome_logs ? '' : "\033[0;95m"
    colorcodes['icyan']      = params.monochrome_logs ? '' : "\033[0;96m"
    colorcodes['iwhite']     = params.monochrome_logs ? '' : "\033[0;97m"

    // Bold High Intensity
    colorcodes['biblack']    = params.monochrome_logs ? '' : "\033[1;90m"
    colorcodes['bired']      = params.monochrome_logs ? '' : "\033[1;91m"
    colorcodes['bigreen']    = params.monochrome_logs ? '' : "\033[1;92m"
    colorcodes['biyellow']   = params.monochrome_logs ? '' : "\033[1;93m"
    colorcodes['biblue']     = params.monochrome_logs ? '' : "\033[1;94m"
    colorcodes['bipurple']   = params.monochrome_logs ? '' : "\033[1;95m"
    colorcodes['bicyan']     = params.monochrome_logs ? '' : "\033[1;96m"
    colorcodes['biwhite']    = params.monochrome_logs ? '' : "\033[1;97m"

    return colorcodes
}

//
// Construct and send completion email
//
def completionEmail(summary_params) {

    // Set up the e-mail variables
    def subject = "[$workflow.manifest.name] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[$workflow.manifest.name] FAILED: $workflow.runName"
    }

    def summary = [:]
    for (group in summary_params.keySet()) {
        summary << summary_params[group]
    }

    def misc_fields = [:]
    misc_fields['Date Started']              = workflow.start
    misc_fields['Date Completed']            = workflow.complete
    misc_fields['Pipeline script file path'] = workflow.scriptFile
    misc_fields['Pipeline script hash ID']   = workflow.scriptId
    if (workflow.repository) misc_fields['Pipeline repository Git URL']    = workflow.repository
    if (workflow.commitId)   misc_fields['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision)   misc_fields['Pipeline Git branch/tag']        = workflow.revision
    misc_fields['Nextflow Version']           = workflow.nextflow.version
    misc_fields['Nextflow Build']             = workflow.nextflow.build
    misc_fields['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    def email_fields = [:]
    email_fields['version']      = NfcoreTemplate.version(workflow)
    email_fields['runName']      = workflow.runName
    email_fields['success']      = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration']     = workflow.duration
    email_fields['exitStatus']   = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport']  = (workflow.errorReport ?: 'None')
    email_fields['commandLine']  = workflow.commandLine
    email_fields['projectDir']   = workflow.projectDir
    email_fields['summary']      = summary << misc_fields

    // Check if we are only sending emails on failure
    def email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine       = new groovy.text.GStringTemplateEngine()
    def tf           = new File("${workflow.projectDir}/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt    = txt_template.toString()

    // Render the HTML template
    def hf            = new File("${workflow.projectDir}/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html    = html_template.toString()

    // Render the sendmail template
    def smail_fields           = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "${workflow.projectDir}" ]
    def sf                     = new File("${workflow.projectDir}/assets/sendmail_template.txt")
    def sendmail_template      = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html          = sendmail_template.toString()

    // Send the HTML e-mail
    Map colors = logColours()
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Sent summary e-mail to $email_address (sendmail)-"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            mail_cmd.execute() << email_html
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Sent summary e-mail to $email_address (mail)-"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }
}

//
// Print pipeline summary on completion
//
def completionSummary() {
    Map colors = logColours()
    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}-"
        } else {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.yellow} Pipeline completed successfully, but with errored process(es) ${colors.reset}-"
        }
    } else {
        log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}-"
    }
}

//
// Construct and send a notification to a web server as JSON e.g. Microsoft Teams and Slack
//
def imNotification(summary_params) {
    def hook_url = params.hook_url

    def summary = [:]
    for (group in summary_params.keySet()) {
        summary << summary_params[group]
    }

    def misc_fields = [:]
    misc_fields['start']                                = workflow.start
    misc_fields['complete']                             = workflow.complete
    misc_fields['scriptfile']                           = workflow.scriptFile
    misc_fields['scriptid']                             = workflow.scriptId
    if (workflow.repository) misc_fields['repository']  = workflow.repository
    if (workflow.commitId)   misc_fields['commitid']    = workflow.commitId
    if (workflow.revision)   misc_fields['revision']    = workflow.revision
    misc_fields['nxf_version']                          = workflow.nextflow.version
    misc_fields['nxf_build']                            = workflow.nextflow.build
    misc_fields['nxf_timestamp']                        = workflow.nextflow.timestamp

    def msg_fields = [:]
    msg_fields['version']      = NfcoreTemplate.version(workflow)
    msg_fields['runName']      = workflow.runName
    msg_fields['success']      = workflow.success
    msg_fields['dateComplete'] = workflow.complete
    msg_fields['duration']     = workflow.duration
    msg_fields['exitStatus']   = workflow.exitStatus
    msg_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    msg_fields['errorReport']  = (workflow.errorReport ?: 'None')
    msg_fields['commandLine']  = workflow.commandLine.replaceFirst(/ +--hook_url +[^ ]+/, "")
    msg_fields['projectDir']   = workflow.projectDir
    msg_fields['summary']      = summary << misc_fields

    // Render the JSON template
    def engine       = new groovy.text.GStringTemplateEngine()
    // Different JSON depending on the service provider
    // Defaults to "Adaptive Cards" (https://adaptivecards.io), except Slack which has its own format
    def json_path     = hook_url.contains("hooks.slack.com") ? "slackreport.json" : "adaptivecard.json"
    def hf            = new File("${workflow.projectDir}/assets/${json_path}")
    def json_template = engine.createTemplate(hf).make(msg_fields)
    def json_message  = json_template.toString()

    // POST
    def post = new URL(hook_url).openConnection();
    post.setRequestMethod("POST")
    post.setDoOutput(true)
    post.setRequestProperty("Content-Type", "application/json")
    post.getOutputStream().write(json_message.getBytes("UTF-8"));
    def postRC = post.getResponseCode();
    if (! postRC.equals(200)) {
        log.warn(post.getErrorStream().getText());
    }
}
