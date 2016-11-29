#Dockerized Tools

Collection of dockerized tools for use on the CGC. 

Tools will be added to the **tools** project on the CGC so that they are easy to find. Apps in the **tools** project can easily be copied to other projects as required.

##Useful commands

Build docker image
```
docker build -t <name> .
```
Run image interactively
```
docker run -ti <name> bash
```
Tag image
```
docker tag <name> cgc-images.sbgenomics.com/wayland/<name>:<tag>
```
Login to CGC repository (use authentication token in place of password)
```
docker login cgc-images.sbgenomics.com
```
Push image to CGC repository
```
docker push cgc-images.sbgenomics.com/wayland/<name>:<tag>
```

##Example 

Assume **Dockerfile** and any other required files are located in a directory called align in our home directory. Change to this directory:
```
cd ~/align
```
Build image according to instructions in **Dockerfile** and give it the name **align**:
```
docker build -t align .
```
Run image interactively to check it is working:
```
docker run -ti align bash
```
Tag image with a version number (v0.0.1) ready for submission to CGC repository:
```
docker tag align cgc-images.sbgenomics.com/wayland/align:v0.0.1
```
Login to CGC repository:
```
docker login cgc-images.sbgenomics.com
```
Push image to CGC repository:
```
docker push cgc-images.sbgenomics.com/wayland/align:v0.0.1
```
