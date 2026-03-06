FROM ubuntu:24.04

# Install basic development tools, ca-certificates, and iptables/ipset, then clean up apt cache to reduce image size
RUN apt-get update && apt-get install -y --no-install-recommends \
  aggregate \
  ca-certificates \
  curl \
  dnsutils \
  fzf \
  nodejs \
  npm \
  gh \
  git \
  gnupg2 \
  iproute2 \
  ipset \
  iptables \
  jq \
  less \
  man-db \
  procps \
  unzip \
  ripgrep \
  zsh \
  tmux \
  vim \
  && rm -rf /var/lib/apt/lists/*

# Install codex
ENV NPM_CONFIG_PREFIX=/usr/local/share/npm-global
RUN mkdir -p $NPM_CONFIG_PREFIX
RUN npm i -g @openai/codex \
  && npm cache clean --force \
  && rm -rf /usr/local/share/npm-global/lib/node_modules/codex-cli/node_modules/.cache \
  && rm -rf /usr/local/share/npm-global/lib/node_modules/codex-cli/tests \
  && rm -rf /usr/local/share/npm-global/lib/node_modules/codex-cli/docs

# Set up non-root user
ARG USERNAME=agent
RUN useradd -m -s /bin/bash $USERNAME
USER $USERNAME
ENV PATH=$PATH:$NPM_CONFIG_PREFIX/bin

# Install Claude Code
RUN curl -fsSL https://claude.ai/install.sh | bash

# Configure git
RUN git config --global user.name "iijimahr" && \
    git config --global user.email "iijimahr@users.noreply.github.com"

