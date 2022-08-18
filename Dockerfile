FROM rust:1.62.1 as builder

COPY . /opt/diploid-contam/
WORKDIR /opt/diploid-contam
RUN cargo install --path .

FROM debian:buster-slim as exec
RUN apt-get update && rm -rf /var/lib/apt/lists/*
COPY --from=build /usr/local/cargo/bin/diploid-contam-estimator /usr/local/bin/diploid-contam-estimator
ENTRYPOINT ["/usr/local/bin/diploid-contam-estimator"]