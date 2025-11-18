function expr_matrix = importfile_tsv(filename, dataLines)
%IMPORTFILE1 从文本文件中导入数据
%  ORIDATA = IMPORTFILE1(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  返回数值数据。
%
%  ORIDATA = IMPORTFILE1(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  oridata = importfile1("D:\Dataset\Source File\Biase\Data\ori_data.tsv", [2, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2023-09-06 10:09:57 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [2, Inf];
end

%% 设置TSV文件路径并读取其列名
% 打开TSV文件进行读取
fileID = fopen(filename, 'r');
% 读取TSV文件的第一行，提取列名
headerLine = fgetl(fileID);
% columnNames = string(split(headerLine));  % tsv:默认; csv:','
columnNames = string(strsplit(headerLine,'\t'));  % 优化 str split
% 关闭文件
fclose(fileID);

% varialbleName:
cell_num = length(columnNames)-1;
columnNames(1) = "Var1";
%% 设置导入选项并导入数据
opts = delimitedTextImportOptions("NumVariables", cell_num+1); % cell num + 1

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% 指定列名称和类型
opts.VariableNames = columnNames;
opts.SelectedVariableNames = columnNames(2:end);
opts.VariableTypes = ["string" repmat("double",1,cell_num)];  %  ['string' repmat('double',1,cell_num)]

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% 导入数据
expr_matrix = readtable(filename, opts);

%% 转换为输出类型
expr_matrix = table2array(expr_matrix);
end