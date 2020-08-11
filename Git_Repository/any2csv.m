function str = any2csv(s,cellsep,humanreadable,file)
%any2text writes contents of s to a nicely formated csv-file
% and/or returns the csv-string. For formating details, see strcells below.
% 
% Inputs
% ------
% s: Any Object, including cells and structs.
%    User-defined objects are tried to be processed as char(object)
% cellsep (optional): separation-characters of cells. Default = ','
% humanreadable (optional): if set to true, empty fields are filled with
%    spaces, so cells are well-aligned in a simple text editor.
% file (optional): file to which csv-string is written. If defined, file is
%    opened with the default application for the specified file(ending).
%
% Outputs
% -------
% str: (1 x n) character-array, newlines = \n. Newlines in csv-file are system-
%    specific.
%
% Example
% -------
% try anything you like, e.g.
%
%    s=struct('a',{1 2 3},'b',rand(4),'c',struct('c1','hello','c2',@(x) 5*x+sin(x)));
%
% for nicely displaying the content, call
%    any2csv(s,'|',1)
%
% for saving the content as csv file:
%    any2csv(s,';',0,'s.csv')
% This command will also launch 's.csv' in the program you specified for
% csv-files in your OS (tested for OS X 10.4.11 and Windows XP SP2)
%
% See also DLMWRITE, CSVWRITE, XLSWRITE.
%
% Version: 1.0, March 2008, Felix Zoergiebel
%          1.1, March 2008, Felix Zoergiebel: empty cells work now
%          1.2, May 2008, Felix Zoergiebel
%               - removed empty structs bug
%               - structs can be aligned vertically or horizontally
%               - content of cell-arrays and struct fields can be
%                 positioned right to and/or under titles/subscripts in
%                 user defined distance.
%               Thanks to Jiro Doke for proposing these features!
%
% Please report any bugs to Felix Zoergiebel (felix_z <at> web.de).
% Feel free to use, modify and improve this script, but please always pass/
% update the version/author history above!

% ----- Init -----
if nargin<2 || isempty(cellsep)
    cellsep=';'; % characters separating cell-contents
end
if nargin<3 || isempty(humanreadable)
    humanreadable=false;
end
cellsep=cellsep(:)';
name=inputname(1);
if isempty(name), name='anonymous variable'; end
% ----------------

% get content of s as cell-array of strings
c=strcells(s);
N=size(c,1);
M=size(c,2);

% add blanks
if humanreadable
    for m=1:M-1 % skip last column
        maxW=max(cellfun(@(x) size(x,2), c(:,m)));
        blanks=repmat(' ',1,maxW);
        for n=1:N
            c{n,m}=[c{n,m} blanks(1:maxW-length(c{n,m}))];
        end
    end
end

% add cell-separators and newlines
for n=1:N
    for m=1:M-1 % skip last column
        c{n,m}=[c{n,m} cellsep];
    end
    c{n,M}=[c{n,M} 10]; % 10 = \n
end

% concatenate all strings in the right order
c=c';
str=['Contents of ' name 10 10 c{:}];

% write to file and open
if nargin==4
    f=fopen(file, 'wt');
    fwrite(f,str,'char');
    fclose(f);
    if isunix
        system(['open ' file]);
    elseif ispc
        system(['start ' file]);
    end
end

end

function cellbox = strcells(s)
%strcells recursively writes contents of s to a cell-array of strings.
% character-arrays: get split into a sinlge cell-column at newlines and at
%    row-endings
% num-/logical-arrays: each value one cell. Arrays with more than 2 dims
%    are split into cell arrays of 2D-arrays. Information about splitting
%    is written to a title-cell
% other array (cell/struct): elements are written to cellboxes in the lay-
%    out of the input-array. Subscript indices are placed above cellboxes.
% struct: fields are displayed with fieldnames as titles.
% everything else: if possible convert to character-array, otherwise
%    display placeholder '*** Object of class %s ***', %s=class(s)

% Layout-adjustments
% ------------------
numformat='%g'; % format string for number conversion
fieldmarker='>> '; % string preceding fieldnames. Size must be 1 x n
              %left  right  top  bottom
gap_c=struct('l',1,'r',0,'t',0,'b',1); % gaps around cell elements
gap_s=struct('l',1,'r',0,'t',0,'b',1); % gaps around struct-fields
gap_n=struct('l',0,'r',0,'t',0,'b',0); % gaps around arrays
gap_t=struct('l',0,'r',0,'t',0,'b',0); % gaps around text

gap_content=struct('l',1,'t',0); % gap between titles and content. Both must be >=0, sum should be >0
vertstruct=false; % align struct vertically or horizontally
% ------------------

siz=size(s);

if ischar(s) || isnumeric(s) || islogical(s)
    if length(siz)<3 % 1D and 2D arrays
        if ischar(s) % string
            cellbox=cell(gap_t.t+gap_t.b+size(s,1),gap_t.l+gap_t.r+1);
            s=[s ones(siz(1),1)*10]'; % append newlines to ends of rows, than flip rows and columns
            eol=[0; find(s==10)]; % find *all* newlines
            for idx=2:length(eol)
                % seperate string at newlines and align parts vertically in cells
                cellbox(gap_t.t+idx-1,gap_t.l+1)={sprintf('"%s"',s(eol(idx-1)+1:eol(idx)-1))};
            end
        elseif isnumeric(s) || islogical(s) % number (array or single value)
            N=siz(1);
            M=siz(2);
            cellbox=cell(N+gap_n.t+gap_n.b,M+gap_n.l+gap_n.r);
            %put each numbers (as string) to a cell
            for n=1:N
                for m=1:M
                    cellbox{n+gap_n.t,m+gap_n.l} = sprintf(numformat,s(n,m));
                end
            end
        end
    else % multidimensional arrays
        % create a cell for all 2-dim arrays in s
        cellsize=siz(3:end);
        subscell=num2cell(char((1:length(siz)-2)+104)); % {i,j,k,...}
        subsstr=sprintf(['%s' repmat(',%s', 1, length(siz)-3)],subscell{:});
        comma1='';
        if length(siz)==3
            cellsize=[cellsize 1];
            comma1=',1';
        end
        multidimcell=cell(cellsize);
        % fill those cells with the 2D arrays
        subs=cell(1,length(siz)-2);
        for ind=1:prod(siz(3:end))
            [subs{:}]=ind2sub(siz(3:end),ind);
            multidimcell{subs{:}}=s(:,:,subs{:});
        end
        % generate infotext about conversion from multidimensional array to cell
        infostring=sprintf('"%gD-array. Displaying A(:,:,%s) as cell(%s%s):"',length(siz),subsstr,subsstr,comma1);
        % get cell-arrays of content
        cellbox=strcells(multidimcell);
        % remove cell-gaps (each 2D-array has its own number-gaps)
        cellbox=cellbox(1+gap_c.t:end-gap_c.b,1+gap_c.l:end-gap_c.r);
        % add infotext
        cellbox=[{infostring} cell(1,size(cellbox,2)-1); cellbox];
    end
elseif prod(siz)>1 % non-(char/num/logical)-arrays (usually cell or struct arrays)
    N=siz(1);
    M=prod(siz(2:end)); % handle all higher dimensions as 2nd dim.
    content=cell(N,M);
    titles=cell(N,M);
    % prepare display of more than two-dimensional subscript indices
    subs=cell(1,length(siz)-1);
    formatstr=['"(%g' repmat(',%g', 1, length(siz)-1) ')"'];
    for n=1:N
        for m=1:M
            % get formated subscripts
            [subs{:}]=ind2sub(siz(2:end),m);
            titles{n,m}={sprintf(formatstr, n, subs{:})};
            % get content
            content{n,m}=strcells(s(n,m));
        end
    end   
    write_titles_and_content(gap_c);
elseif isstruct(s) % single struct
    % process similar to cell-array, just replace subscripts by fieldnames
    fields=fieldnames(s);
    L=length(fields);
    if L>0
        if vertstruct, M=1; N=L;
        else M=L; N=1; end
        content=cell(N,M);
        titles=cell(N,M);
        % get content and fieldnames
        for n=1:L
            titles{n}={[fieldmarker fields{n}]};
            content{n}=strcells({s.(fields{n})});
        end
        write_titles_and_content(gap_s);
    else
        cellbox=cell(gap_s.t+1+gap_s.b,gap_s.l+1+gap_s.r);
        cellbox(gap_s.t+1,gap_s.b+1)={'"struct with no fields"'};
    end
elseif iscell(s) % single cells
    if ~isempty(s)
        cellbox=strcells(s{:}); % simply pass content
    else
        cellbox={'"empty cell"'};
    end
else % everything else, e.g. function handles
    try
        str = char(s);
    catch
        str = sprintf('*** Object of class %s ***', class(s) );
    end
    cellbox = strcells(str);
end

% subfunction for cells-arrays and structs
function write_titles_and_content(gap_x)
    % maximum width of column m plus gap right
    maxWc=arrayfun(@(m) max(cellfun(@(x) max(size(x,2),gap_content.l), content(:,m))) + gap_x.r, 1:M);
    % maximum height of row n plus gap bottom
    maxHc=arrayfun(@(n) max(cellfun(@(x) max(size(x,1),gap_content.t), content(n,:))) + gap_x.b, 1:N);
    Ht=gap_content.t+gap_x.t;
    Wt=gap_content.l+gap_x.l;
    cellbox=cell(sum(maxHc+Ht), sum(maxWc+Wt));
    for b=1:N
        for a=1:M
            % write titles
            nt=sum(maxHc(1:b-1)+Ht)+gap_x.t+1;
            mt=sum(maxWc(1:a-1)+Wt)+gap_x.l+1;
            cellbox(nt,mt)=titles{b,a};
            % write contents
            [hc,wc]=size(content{b,a});
            nac=nt+gap_content.t;
            nbc=nac+hc-1;
            mac=mt+gap_content.l;
            mbc=mac+wc-1;
            cellbox(nac:nbc,mac:mbc)=content{b,a};
        end
    end
end


end