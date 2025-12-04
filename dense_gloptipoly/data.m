function data_tens_dict = data(path)
    %current_file = mfilename('fullpath');   % obtain the path of current files
    %[current_path, ~, ~] = fileparts(current_file);   %obtain the path of current directory
    %[parent_path, ~, ~] = fileparts(current_path);  %obtain the path of the parent directory
    
    %data_path = "..\assets\data\literature\";
    data_path = path;
    data_files = dir(fullfile(data_path,'*.csv'));
    data_names = {data_files.name}';
    
    tens_data_cp = cell(length(data_names),1);
    for i=1:length(data_names)
        fulltrace = strcat(data_path,data_names{i});
        A = importdata(fulltrace);% A.data: nonzero values of tensor; 
                                  % A.textdata: subscripts of nonzero
                                  % values

        A_name = data_names{i}; % 'ex02_n10_d3_17Fan-Exa4.3_.csv'
        spl = split(A_name,'_');
        n = eval(spl{2}(2:end));% dimension of tensor
        d = eval(spl{3}(2:end));% order of tensor

        subscripts = get_subscripts(d,n);
        observes = zeros(size(subscripts,1),1);

        if length(A.data) == size(subscripts,1)
            observes = A.data;
        else
            for j = 1:length(A.data)
                index = findrow(eval(A.textdata{j+1}),subscripts);
                observes(index) = A.data(j);
            end
        end
        tens_data_cp{i}.subscripts = subscripts;
        tens_data_cp{i}.observes = observes;
    end
    
    data_tens_dict = containers.Map(data_names,tens_data_cp);
end