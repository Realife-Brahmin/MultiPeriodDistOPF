% Is MA57 usable through MATLAB on this machine, and does MATLAB really use it?
%
% MATLAB ships HSL's MA57 as bin/<arch>/libmwma57.dll under MathWorks' licence and
% uses it for sparse symmetric indefinite LDL' (ldl, decomposition(...,'ldl') and
% backslash). This probe checks the licence, the library, and -- via spparms'
% diagnostic output -- that backslash on a small KKT matrix reports MA57.
%
%   matlab -batch "run('ddp/examples/power_system/ma57_availability_probe_matlab.m')"

out = fullfile('ddp', 'results', 'kkt_ordering', 'ma57_availability_matlab.txt');
if exist(out, 'file'), delete(out); end
diary(out);
fprintf('MATLAB %s on %s\n', version, computer('arch'));
fprintf('license(''test'',''MATLAB'') = %d\n', license('test', 'MATLAB'));
fprintf('maxNumCompThreads = %d\n', maxNumCompThreads);
lib = fullfile(matlabroot, 'bin', computer('arch'), 'libmwma57.dll');
fprintf('libmwma57.dll present: %d  (%s)\n', exist(lib, 'file') == 2, lib);

rng(1);
n = 300; m = 120;
H = sprandsym(n, 0.02) + speye(n);
A = [speye(m), sprand(m, n - m, 0.05)];
K = [H, A'; A, sparse(m, m)];
b = ones(n + m, 1);

fprintf('\n--- backslash with spparms(''spumoni'',2) ---\n');
spparms('spumoni', 2);
x = K \ b;
spparms('spumoni', 0);
fprintf('backslash relative residual %.3e\n', norm(K * x - b) / norm(b));

fprintf('\n--- ldl and decomposition ---\n');
[L, D, P, S] = ldl(K, 'vector');
fprintf('ldl: nnz(L) = %d, 2x2 blocks in D = %d\n', nnz(L), nnz(diag(D, 1)));
dK = decomposition(K, 'ldl');
disp(dK);
fprintf('decomposition relative residual %.3e\n', norm(K * (dK \ b) - b) / norm(b));
diary off;
