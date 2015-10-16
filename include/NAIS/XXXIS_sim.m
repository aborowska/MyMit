function theta_sim = IS_sim(S,n,theta_smooth,par_KFS,RND)
% sample S/2 trajectories for the signal using e.g. the simulation smoother
    S = S/2;
    theta_sim = zeros(n,S);
    for ii = 1:S
        theta_sim(:,ii) = SimSmooth(theta_smooth,par_KFS,RND(:,ii*[2,2]-[1,0]));
    end
end