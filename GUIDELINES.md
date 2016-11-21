# Good coding guidelines
https://tedvinke.wordpress.com/2015/03/15/basic-groovy-and-grails-code-review-guidelines/

## General
- Avoid abbreviations.

## Process

### Naming convention
The name should begin with a capital letter and be descriptive of what the Process does. A verb in present tense followed by a noun is prefered. The noun should be singular if a single process is executed for each element (parrallelization). The noun should be plural if a single process is executed for several elements. Each different elements of the name should begin with a capital letter.
example:
-`DoThing`: this process does one thing
-`DoThings`: this process does multiple things (could be the same thing but multiple times)
If the process is used to execute a specific tool, its name should follow the naming convention.
For example, with the tool `triDent`:
- `ExecuteTrident`

### tag
Use an appropriate tag to easylly follow wich occurence of the process is currently running.
The tag should be unique for each occurence of the process.

## Channel
