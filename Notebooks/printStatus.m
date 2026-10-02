[okC, commit] = system('git -C "." rev-parse --short HEAD');
[okB, branch] = system('git -C "." rev-parse --abbrev-ref HEAD');
[~,   dirty]  = system('git -C "." status --porcelain --untracked-files=no');
if okC == 0 && okB == 0
    if isempty(strtrim(dirty)), state = ''; else, state = ', working tree modified'; end
    fprintf('Run for commit %s on branch %s%s\n', strtrim(commit), strtrim(branch), state);
else
    fprintf('Run from %s (git commit unavailable)\n', pwd);
end
fprintf('Run on %s\n', string(datetime('now','Format','yyyy-MM-dd HH:mm')));