FROM julia:1.1.0

ENV HOME /home/developer
WORKDIR $HOME

USER developer

RUN julia -e "using Pkg; Pkg.add(\"HTTP\");Pkg.add(\"JSON\")"

COPY runModel.ipynb $HOME
COPY InputData.xlsx $HOME