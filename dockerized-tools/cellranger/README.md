# Instructions for building cellranger docker image

If cell profiler is added in one step, then the layer created (1.8GB) is too large to upload to the CGC repository. Therefore we have to add the software in stages to several different layers.

## Step 1: get the required files
Put Dockerfile in an empty directory and download cellranger tarball to this directory:
```
curl -ko cellranger-1.3.1.tar.gz "https://s3-us-west-2.amazonaws.com/10x.downloads/cellranger-1.3.1.tar.gz?AWSAccessKeyId=AKIAJAZONYDS6QUPQVBA&Expires=1488027285&Signature=ECmZCuIwtpns0MeMvndnZFEhMCg%3D"
```

unpack cellranger tarball
```
tar -xzf cellranger-1.3.1.tar.gz
```

## Step 2: break-down cellranger into components that are within the size limit of an image layer
move the larger components of cellranger out of the cellranger-1.3.1 directory:
```
mv cellranger-1.3.1/anaconda-cr-cs .
mv cellranger-1.3.1/cellranger-cs .
mv cellranger-1.3.1/cellranger-tiny-ref .
```

## Step 3: build docker image
```
docker build -t cellranger .
```

## Step 4: tag docker image
```
docker tag cellranger cgc-images.sbgenomics.com/wayland/cellranger:v0.0.3
```

## Step 5: login in to repository
```
login cgc-images.sbgenomics.com
```

## Step 6: push image to repository
```
docker push cgc-images.sbgenomics.com/shapni/cellranger:v0.0.3
```

