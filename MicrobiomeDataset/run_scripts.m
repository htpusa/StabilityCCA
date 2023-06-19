if ~isfolder('data')
    download_data;
end
if ~isfile('data.mat')
    process_data;
end

mkdir('res')
runPipe("species","idMetabolites")
%runPipe("genus","idMetabolites")
runPipe("enzymes","idMetabolites")
runPipe("species","metabolites")
%runPipe("genus","metabolites")
runPipe("enzymes","metabolites")

if ~isfolder('plots')
    do_plots;
end

if ~isfolder('tables')
    collect_results;
end