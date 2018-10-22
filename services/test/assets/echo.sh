#!/bin/bash

sleep 10
echo -e "STDERR output\t with Tab" > /dev/stderr
sleep 5
echo -e "STDERR output\nWith newline" > /dev/stderr
cat /example.sam > /dev/stdout
