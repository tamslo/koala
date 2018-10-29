#!/bin/bash

echo -e "STDERR output\t with Tab" > /dev/stderr
sleep 2
echo -e "STDERR output\nWith newline" > /dev/stderr
cat /example.sam > /dev/stdout
