FROM ubuntu:16.04

ENV NOVO_FILE novocraftV3.09.00.Linux3.10.0.tar.gz
ENV NOVO_LICENSE novoalign.lic

ENV NOVO_PATH /opt/novocraft
RUN mkdir $NOVO_PATH
WORKDIR $NOVO_PATH

# Copy Novocraft suite and NovoAlign license if present
COPY ./assets/$NOVO_FILE ./assets/$NOVO_LICENSE* $NOVO_PATH/

# Extract Novocraft suite to $NOVO_PATH
RUN tar -xzf $NOVO_FILE -C /opt
WORKDIR /

# Add to PATH
ENV PATH $NOVO_PATH:$PATH
