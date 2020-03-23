FROM julia:1.1.1

ENV HOME /home
WORKDIR $HOME

COPY runModel.ipynb $HOME
COPY InputData.xlsx $HOME


