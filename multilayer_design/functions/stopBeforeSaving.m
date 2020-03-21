function stopBeforeSaving(name)
% the warn is actually used for waiting matlab to plot
w=warndlg("Check Command line instruction","warn");
close(w)

input(strcat("\nDo you wand to save the file as '",name,...
                 "' ? \n Press enter to continue",...
                 "\n otherwize stop the excecution otherwise."));
    
end