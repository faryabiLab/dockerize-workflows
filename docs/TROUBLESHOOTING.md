# Common Errors
Here are some common errors one might run into when running a workflow, and steps to solve them.

## Failed to Connect to Database 
Example error message:
```
[2025-10-16 11:26:33,35] [info] Running with database db.url = jdbc:mysql://0.0.0.0:52000/PipelineDatabase?allowPublicKeyRetrieval=true&rewriteBatchedStatements=true&useSSL=false
[2025-10-16 11:26:43,54] [error] Failed to instantiate Cromwell System. Shutting down Cromwell.
java.sql.SQLTransientConnectionException: db - Connection is not available, request timed out after 10000ms.
```
Solution: Ensure the MySQL database Docker image is currently running. 
1. Run `docker ps` on the same machine you executed the workflow from.
2. The MySQL image entry should look something like:
```
80420dd85e5a   mysql/mysql-server:latest   "/entrypoint.sh mysq…"   12 hours ago        Up 12 hours (healthy)   33060-33061/tcp, 0.0.0.0:52000->3306/tcp, [::]:52000->3306/tcp   PipelineDatabase
```
3. If the database cannot conenct, this means you should not see this line in the output. Run `workflows/start_database.sh` to start the MySQL Docker image.
4. Check periodically with `docker ps` until the `STATUS` column turns from `Up <uptime> (starting)` to `Up <uptime> (healthy)`.
5. Restart your workflow.

## Java Runtime Version
Example error message:
```
Uncaught error from thread [cromwell-system-akka.dispatchers.engine-dispatcher-6]: wdl/draft3/parser/WdlParser$Ast has been compiled by a more recent version of the Java Runtime (class file version 61.0), this version of the Java Runtime only recognizes class file versions up to 55.0, shutting down JVM since 'akka.jvm-exit-on-fatal-error' is enabled for ActorSystem[cromwell-system]
...
```
In short, this is caused by running Cromwell with too old a Java version. In the above example, we tried to run the Cromwell engine jarfile with Java 8 (class file version 55.0, Java 8), but this version of Cromwell requires class file version 61.0 (Java 17). 

Most likely, the root cause is having a conda environment loaded with an old Java version while running the workflow. If this is the case, simply deactivate with `conda deactivate` and try again. 

If not such environment is activated (and you cannot activate a conda environment that does have the correct Java version), update your Java with `sudo apt-get update
&& sudo apt-get install openjdk-17-jre`.

If future version conflicts occur, mappings between class file versions and Java versions can be found [here](https://javaalmanac.io/bytecode/versions/).

## Logging on an NFS
This isn't so much a fatal error as it is a point of confusion. When running cromwell manually (i.e. `cromwell.jar run`), it's often desirable to redirect output to a file such as `nohup.out` instead of having `stdout`revieve the output, which relies on the terminal session staying alive.

Logging on an NFS is fine, but due to automatic file handle recycling, this file may become unmounted during longer jobs. This will _not_ break the workflow, but workflow logging will be directed to a `.nfsXXXXXX` file instead of `nohup.out`. The workflow will continue executing as normal.
