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

# Install required packages
RUN apt-get update && \
    apt-get install -y \
    ipython \
    python-lxml \
    python-matplotlib \
    python-nose \
    python-pip* \
    python-scipy 

RUN pip install --upgrade pip
RUN pip install --user --no-cache-dir notebook==5.*
RUN pip install --user ipykernel==4.7.0

USER developer

ENV HOME /home/developer
WORKDIR $HOME

ENV PATH $PATH:$HOME
ENV PATH $PATH:$HOME/.local/bin

COPY runModel.ipynb $HOME
COPY InputData.xlsx $HOME


