# Cloud Computing

## Objectives
* What is it?
* What are the tradeoffs of cloud computing?
* What are the benefits?

## Verifying your environment:

```
whoami # Provides username

df -h # Space on hard drive

cat /proc/cpuinfo # Details on how many processors the machine has

tree -L 1 # Shows diagram of the file system 1 level below your current location

```

## Staying connected:
screen and tmux are programs that keep jobs running even when you are disconnected

```
tmux new -s descriptive_name

control b d # to detatch

tmux list-sessions

tmux attach -t descriptive_name # to connect to an existing session

tmux kill-session -t session_name # kill a session

```

## Install additional software:
* tmux is not automatically installed
* Project Managers: YUM (Red Hat and Centos instances) or APT (Debian or Ubuntu instances)

```
apt search tmux # search for packages

# Try using it to search for your favorite bioinformatics program:
apt search blast


# We can't do this because we don't have administrative privileages

sudo apt upgrade # run apt from the adminsitrator account
apt install PackageName
```

