FROM rust:bullseye AS builder

RUN apt-get -y update \
    && apt-get -y dist-upgrade \
    && apt-get install -y --no-install-recommends --no-install-suggests libclang-dev

WORKDIR /usr/src

RUN cargo new --bin gc_prop_in_window
WORKDIR /usr/src/gc_prop_in_window

COPY Cargo.toml Cargo.lock ./
RUN cargo build --release
RUN rm -rf src

COPY src ./src
RUN rm ./target/release/deps/gc_prop_in_window*
RUN cargo build --release

FROM rust:bullseye

WORKDIR /app

COPY --from=builder /usr/src/gc_prop_in_window/target/release/gc_prop_in_window .
