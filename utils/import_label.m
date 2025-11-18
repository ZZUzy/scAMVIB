function label = import_label(filename, dataLines)
%IMPORTFILE1 从文本文件中导入数据
%  SUBTYPE = IMPORTFILE1(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  返回数值数据。
%
%  SUBTYPE = IMPORTFILE1(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  subtype = importfile1("D:\Experiment\Code\Python Code\scGAC-main\data\Buettner\label.tsv", [2, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2023-09-06 10:32:14 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [2, Inf];
end

%% 设置导入选项并导入数据
opts = delimitedTextImportOptions("NumVariables", 2);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% 指定列名称和类型
opts.VariableNames = ["Var1", "VarName2"];
opts.SelectedVariableNames = "VarName2";
opts.VariableTypes = ["string", "string"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% 指定变量属性
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% 导入数据
label = readtable(filename, opts);

%% 转换为输出类型
label = table2array(label);
end