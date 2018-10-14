#!/bin/bash

sleep 10
echo -e "STDOUT output\t with Tab" > /dev/stdout
echo -e "STDERR output\t with Tab" > /dev/stderr
sleep 5
echo -e "STDOUT output\nWith newline" > /dev/stdout
echo -e "STDERR output\nWith newline" > /dev/stderr
