# Pull and run docker image

The image contains R studio server, R 3.6.1, bioC 3.10, Mitch 0.2.3, DESeq2, 
edgeR, limma, ABSSeq, Seurat, Sleuth, Fishpond(swish), Muscat

```
docker pull mziemann/mitch
docker run -it -v ${PWD}:/mnt -p 8787:8787 mziemann/mitch
```

Then you have the choice to run R from the command line or connect to the 
R-studio server.

On the docker command line run:
```
service rstudio-server start
```

Credentials are mitch:mitch

To ensure that R studio remains avalable, the docker run command can be 
executed in a tmux shell.

Connect to Rstudio by pointing the webbrowser to the IP and specify the port as 
8787. For example if running on a local machine use localhost:8787
