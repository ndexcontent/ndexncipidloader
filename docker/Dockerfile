# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get install -y build-essential 

RUN mkdir /tmp/ndexncipidloader
COPY ./ /tmp/ndexncipidloader/
RUN pip install /tmp/ndexncipidloader

RUN rm -rf /tmp/ndexncipidloader

RUN mkdir /ndexncipidloader
WORKDIR /ndexncipidloader

COPY docker/README.txt /ndexncipidloader/

CMD /bin/cat /ndexncipidloader/README.txt
