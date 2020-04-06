clear, close all

files=split(string(ls('designs/*.mat')))
for i=1:length(files)
    
    design_file = files(i)
    load(design_file)
    d_layers(1)=1e-6;
    d_layers(end)=3e-6;
    save(design_file,'idx_layers','d_layers','n_eff1','n_eff2','n_eff3')
end