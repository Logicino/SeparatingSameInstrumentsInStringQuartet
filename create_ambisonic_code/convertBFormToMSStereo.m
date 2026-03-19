function stereoOut = convertBFormToMSStereo(WXYZ)
% B-format IR (WXYZ) and convert to mid-side stereo (stereoOut)

%formulation
%Zotter, F. and Frank, M., 2019. Ambisonics: A practical 3D audio theory for recording, studio production, sound reinforcement, and virtual reality (p. 210). Springer Nature.
% Eqn 1.7

alpha = 0.5;

W = WXYZ(:,1);
Y = WXYZ(:,3);

M = 0.5*[1,1;1,-1]*[2-alpha,0;0,alpha];

stereoOut = zeros(length(W),2);
for n=1:length(W)
    stereoOut(n,:) = M*[W(n);Y(n)];
end

end 