FROM julia:1.1.1

# Replace 1000 with your user / group id
RUN export uid=1000 gid=1000 && \
    mkdir -p /home/developer && \
    mkdir -p /etc/sudoers.d && \
    echo "developer:x:${uid}:${gid}:Developer,,,:/home/developer:/bin/bash" >> /etc/passwd && \
    echo "developer:x:${uid}:" >> /etc/group && \
    echo "developer ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/developer && \
    chmod 0440 /etc/sudoers.d/developer && \
    chown ${uid}:${gid} -R /home/developer

RUN apt-get update && \
	apt-get install python-pip -y

USER developer

RUN pip install --user --no-cache-dir notebook==5.*
RUN pip install --user ipykernel==4.7.0

ENV HOME /home/developer
WORKDIR $HOME

ENV PATH $PATH:$HOME
ENV PATH $PATH:$HOME/.local/bin

RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add Dates     	;precompile");using Dates'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add DataFrames	;precompile");using DataFrames'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add XLSX       	;precompile");using XLSX'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add CSV        	;precompile");using CSV'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add Missings   	;precompile");using Missings'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add ControlSystems;precompile");using ControlSystems'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add MAT      		;precompile");using MAT'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add Statistics    ;precompile");using Statistics'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add LinearAlgebra ;precompile");using LinearAlgebra'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add Distributions ;precompile");using Distributions'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add JuMP    		;precompile");using JuMP'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add Ipopt	    	;precompile");using Ipopt'
RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add Plots	    	;precompile");using Plots'

COPY Integrated_model.jl $HOME
COPY InputData.xlsx $HOME