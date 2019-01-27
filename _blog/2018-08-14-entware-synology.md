---
layout: posts
title: "Setting up and using Entware on Synology device"
date: 2018-08-14
excerpt: "Use Entware package manager on a Synology storage device."
---

Taken from [https://keestalkstech.com/2018/03/install-nano-with-entware-on-synology-nas-dsm6/](https://keestalkstech.com/2018/03/install-nano-with-entware-on-synology-nas-dsm6/).

1. Prepare directory for Entware

```bash
# create a home for Entware
mkdir -p /volume1/@Entware/opt
# go on as root
sudo -i
# remove a previous install
rm -rf /opt
# link the folders
ln -sf /volume1/@Entware/opt /opt
echo "Done!"
```

2. Determine system architecture

```bash
\
printf "\nProcessor:   "; \
cat /proc/cpuinfo | \
grep "model name" | \
grep "[^:]*$" -o  | \
uniq; \
printf "Architecture: "; \
uname -m; \
printf "\n"
```

3. Install Entware by script - use correct script based on architecture

```bash
# armv5
wget -O - http://bin.entware.net/armv5sf-k3.2/installer/generic.sh | /bin/sh
# armv7
wget -O - http://bin.entware.net/armv7sf-k3.2/installer/generic.sh | /bin/sh
# armv8
wget -O - http://bin.entware.net/aarch64-k3.10/installer/generic.sh | /bin/sh
# x64
wget -O - http://bin.entware.net/x64-k3.2/installer/generic.sh | /bin/sh
```

4. Enable Entware to startup by default (DS6+ installation)

```bash
# leave root
exit;
# remove previous file
rm entware-startup.sh 2> /dev/null
# write the startup file
printf "#!" >> entware-startup.sh
echo "/bin/sh" >> entware-startup.sh
echo "" >> entware-startup.sh
echo "case $1 in" >> entware-startup.sh
echo "    start)" >> entware-startup.sh
echo "    mkdir -p /opt" >> entware-startup.sh
echo "    mount -o bind /volume1/@Entware/opt /opt" >> entware-startup.sh
echo "    /opt/etc/init.d/rc.unslung start" >> entware-startup.sh
echo "    ;;" >> entware-startup.sh
echo "    stop)" >> entware-startup.sh
echo "    ;;" >> entware-startup.sh
echo "esac" >> entware-startup.sh
# copy the startup file
sudo mv entware-startup.sh /usr/local/etc/rc.d/entware-startup.sh
echo "Done!"
```

5. Enable script and reboot system

```bash
sudo -i
echo "" >> /etc/profile;
echo ". /opt/etc/profile" >> /etc/profile
reboot
```

6. Install Entware packages

```bash
# example with nano (modify accordingly for other packages)
sudo opkg install nano
```

Here are all of the packages available for x86-64 architectures:

[https://pkg.entware.net/binaries/x86-64/Packages.html](https://pkg.entware.net/binaries/x86-64/Packages.html)

To view packages for other architectures, visit following link and navigate to architecture and `Packages.html` file:

[https://pkg.entware.net/binaries/](https://pkg.entware.net/binaries/)
