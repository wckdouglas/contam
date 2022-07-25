FROM rust:1.62.1 as builder

COPY . /opt/diploid-contam/
WORKDIR /opt/diploid-contam
RUN cargo install --path .

FROM builder as final
ENV RUST_LOG=info
ENTRYPOINT ["/usr/local/cargo/bin/diploid-contam-estimator"]
