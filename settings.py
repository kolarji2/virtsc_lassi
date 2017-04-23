#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging
local_weight_function=''
global_weight_function=''
morgan_radius=2
threshold=1e-9

def init_loggging():
	logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
