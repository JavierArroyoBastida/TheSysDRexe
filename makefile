IMG_NAME=jl

COMMAND_RUN=docker run \
	  --name ${IMG_NAME} \
	  --detach=true \
	  --rm \
 	  -it \
	  ${IMG_NAME} 

build_jl_image:
	docker build --no-cache --rm -t ${IMG_NAME} .

remove_jl_image:
	docker rmi ${IMG_NAME}

run_jl:
	$(COMMAND_RUN)

stop_jl:
	docker stop ${IMG_NAME}

	
