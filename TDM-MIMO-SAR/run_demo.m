% 一键运行入口：在 MATLAB 中打开本文件，点击“运行”。
% 所有配置集中在 default_config.m；输出写入 results 子目录。
clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));
out = main_tdm_mimo_sar;
report = validate_simulation(out);
