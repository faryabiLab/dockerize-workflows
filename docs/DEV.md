## Development Guide
Notes and tips on future development of the workflows for other developers.

### GitHub Actions
This repository already has several GitHub Actions set up. These are:
1. Auto - building of `imports.zip` from `wdl_tasks/`
2. Auto - checking the validity of the WDL workflowsd woth `womtool`
3. Building useful graphical representations of each workflow

As such, when developing locally, remember to:
* Modify task files in `wdl_tasks/`, push them, and then `git pull` after a few seconds to retrieve the updated `imports.zip` in each workflow directory. Much faster than manually zipping the directory every time; and makes sure that changes persist on GitHub.
* Check the Actions tab for the workflow validation check; this will quickly tell you whether or not the workflow is valid in the eyes of Cromwell.
