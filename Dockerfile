FROM rust:1.62.1 as builder

COPY . /opt/diploid-contam/
WORKDIR /opt/diploid-contam
RUN cargo install --path .

FROM builder as final
CMD ["diploid-contam-estimator"]