function  [transmitMat,compParam] = lz77(inputVar,windowType,windowSize)
%Author :           Vishnu Muralidharan A25208488
%                   Dept. of Electrical Engineering, 
%                   Univ. of Alabama in Hutsville
%                   Done for EE 610:Digital Signal Compression
%                   A Sliding Window Based Algorithm for Compression of Data 
%****************************************************************************
%Inputs:            inputVar = input variable to be compressed. can be a string or a text file
%                   windowType = decides the window type for the dictionary
%                   used for compression: 'Fixed' or 'Growing'
%                   windowSize = number that specifies size of the sliding
%                   window
%Outputs:           transmitMat = table variable that returns the compression result 
%                   compParam = struct variable that returns the
%                   compression paramaters
if ismember('.txt',inputVar)
    str = [inputVar];
    formatSpec = '%s';
    fileID = fopen(str,'rt');
    input = fscanf(fileID,formatSpec);
else
    input = inputVar;
end
    
switch nargin
    case 1
        dictType = 'Growing';
    
    case 2
        dictType = windowType;
        if strcmp(dictType,'Fixed')
            error('Window Size is not specified');
        end
        
    case 3
        dictType = windowType;
        if strcmp(dictType,'Fixed')
            dictSize = windowSize;
        end
        
end
            
%% initialization of variables

len = 0;
position = [];                  % vector that stores the psoition elements
matchLen = [];                  % vector that stores the length of match
lbPtr = 1;                      % lookahead buffer pointer
sBuf = [];                      % seacrh buffer
symbol = [];                    % vector that stores symbol elements to transmit

%% Begin compression process

while ~isempty(input(:,lbPtr:end))          % loop till lookahead buffer is not empty
    if isempty(sBuf)                        % for intial iteration
        position = [position;0];            % position = 0 since search buffer pointer is uninitialized
        matchLen = [matchLen;0];            % matched length = 0 since search buffer pointer is uninitialized
        symbol = [symbol;input(lbPtr)];     % symbol = current element since search buffer pointer is uninitialized
        sbPtr = 1;                          % initialize seacrh buffer pointer to first element in imput
        lbPtr = lbPtr + 1;                  % increment lookahead buffer pointer
        sBuf = input(:,1:sbPtr);            % store array till search buffer pointer
    else
        currSym = input(lbPtr);
        idxSym = lbPtr;                     % get current symbol to be processed
        idx = find(sBuf == input(lbPtr));   % find matching symbols for current symbol in search buffer
        if numel(idx) ~= 0

            posTemp = idxSym - idx;                     %get offsets of matching symbols
            lenArr = zeros(1,length(idx));              %array that stores match lengths for macthed symbols
            matchPos = 1;
            sbTemp = idx(matchPos);                     %intialize temporary pointer to first symbol match position
            lbTemp = idxSym;                            %initialize temporary pointer to current symbol to encoded 
            while 1
                % for each matching symbol in search buffer, find symbol
                % match length
                if input(sbTemp) == input(lbTemp)
                    %if match found, increment length of match, sbtemp and
                    %lbTemp to check for next elements
                    lenArr(matchPos) = lenArr(matchPos) + 1;
                    sbTemp = sbTemp + 1;
                    lbTemp = lbTemp + 1;
                    if lbTemp > size(input,2)
                        break;
                    end
                else
                    matchPos = matchPos + 1;
                    if matchPos > length(idx)
                        break;
                    end
                    sbTemp = idx(matchPos);
                    lbTemp = idxSym;
                end
            end
            maxLen = max(lenArr);                           %find largest symbol match length
            matchIdx = max(find(lenArr == max(lenArr)));    %find symbol position in search buffer for longest match
            pos = posTemp(matchIdx);                        %get offset of longest symbol match
            symbol = [symbol;'O'];                          %append null symbol since match was found
        else
            % if no match found, set position to 0, length to 0 and
            % transmit current symbol
%             fprintf('a matching symbol was not found\n');
            pos = 0;                                        
            maxLen = 0;
            symbol = [symbol;input(idxSym)];
        end
        position = [position;pos];
        matchLen = [matchLen;maxLen];
        %if symbol was matched to L bytes, next symbol to be encoded is L
        %bytes forward
        %else next symbol in lookahead buffer needs to be processed
        if maxLen == 0
            lbPtr = lbPtr + 1;
        else
            lbPtr = idxSym + maxLen;
        end
        %set sbPtr to lbPtr - 1 since the search buffer is till last element pointed by lbPtr 
        sbPtr = lbPtr - 1;
        if exist('dictSize','var')
            if dictSize > sbPtr - dictSize
                sBuf = input(:,1:sbPtr);                        % update search buffer
            else
                sBuf = input(:,sbPtr-6:sbPtr);                        % update search buffer
            end
        else
            sBuf = input(:,1:sbPtr);                        % update search buffer
        end
    end
end
transmitMat = table(position,matchLen,symbol);
compOut = table2struct(transmitMat);
%% Calculating bits required and compression ratio

% Number of bits to convert position array to binary
sumbit = 0;
for i=1:1:length(position)
    sumbit = sumbit + length(de2bi(position(i)));
end

% Number of bits to convert length array to binary
sumbitlen = 0;
for i=1:1:length(matchLen)
    sumbitlen = sumbitlen + length(de2bi(matchLen(i)));
end

% the codewords for the symbols are egenarted using Huffman code

transSymbols = unique(symbol);          % get source alphabet

%get probabilities for each symbol in alphabet
for i=1:1:length(transSymbols)
    occur = find(symbol == transSymbols(i));
    prob(i) = length(occur)/length(symbol);
end

% convert to cell to use huffmandict function
transSymbols = cellstr(unique(symbol));
[dict,avglen] = huffmandict(transSymbols,prob);
codewords = huffmanenco(cellstr(symbol),dict);
orgBits = size(input,2)*8;
compBits = sumbit+sumbitlen+length(codewords);
fprintf('Total number of bits before compression:\t\t %d bits\n',orgBits);
fprintf('Total number of bits after compression:\t\t\t %d bits\n',compBits);
orgBits = size(input,2)*8;
compRatio =orgBits/compBits;
fprintf('Compression ratio:\t\t\t\t\t\t\t\t %f\n',compRatio);
compParam.originalBits = orgBits;
compParam.compressedBits = compBits;
compParam.compressionRatio = compRatio;

%% Decoder
output =[];                         % initialize decompressed bitstream
for i=1:1:size(position,1)
    if symbol(i) ~= 'O'                 % if a non-null symbol was received
        output = [output;symbol(i)];    % append symbol to output
    else
        offset = size(output,1) - position(i) + 1;          %if a null symbol was received
        if offset+matchLen(i)-1 < size(output,1)
            appendArr = output(offset:offset+matchLen(i)-1);    % go back in output to position specified by offset, count symbols equal to length and append array to output
            output = [output;appendArr];                        %apopend decoded array to deocded bitstream
        else
            excess = offset+matchLen(i)-1 - size(output,1);     %in certain cases, symbols form lookahead buffer are to be copied
            appendTemp = output(offset:end);                    % refer to book for decoding <3,5,'O'>
            outputTemp = [output;appendTemp];
            offsetTemp = size(outputTemp,1)-excess;
            appendArr = outputTemp(offsetTemp:offsetTemp+excess-1);
            output = [outputTemp;appendArr];
        end
    end
end
output = output';

%% Lossless check
if isequal(output,input)
    fprintf('the sequence was decoded correctly\n');
    compParam.losslessChk = 1;
else
    fprintf('the sequence was not decoded correctly\n');
    compParam.losslessChk = 0;
end