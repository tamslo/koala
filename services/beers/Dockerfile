# Downloads beers_simulator, execute with `reads_simulator.pl [options]`

FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y git g++ make libswitch-perl ruby-full

RUN cpan Math::Random

# Install simulator
ENV DST=/opt/beers_simulator
RUN git clone https://github.com/tamslo/beers_simulator.git $DST
COPY ./assets/make_config_files_for_subset_of_gene_ids.pl $DST
RUN chmod a+x $DST/reads_simulator.pl
ENV PATH $DST:$PATH

# Install evaluation
RUN gem install erubis
ENV DST=/opt/beers_evaluation
RUN git clone https://github.com/khayer/aligner_benchmark.git $DST
ENV PATH $DST:$PATH
