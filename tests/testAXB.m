% we generate test matrices with EIDORS http://eidors.org in Matlab
imdl=mk_common_model('a2C',4);
img=mk_image(imdl,1);
img.elem_data = img.elem_data - elem_select(img.fwd_model,'x>0 & y<0.25')/2;
img.elem_data = img.elem_data + elem_select(img.fwd_model,'x<0.25 & y>0')/3;
clf; show_fem(img,1);
img.fwd_solve.get_all_nodes= 1;
A=system_mat_1st_order(img); A=A.E;
A(1,:)=0; A(:,1)=0; A(1,1)=1; % node1 is gnd node
B=[img.fwd_model.stimulation(:).stim_pattern];
B=[zeros(size(img.fwd_model.nodes,1),length(img.fwd_model.stimulation));B];
X=fwd_solve(img); X=X.volt;
Xx = A\B;
tol=max(max(abs(X-Xx)));
fprintf('A\\B == fwd_solve(img) ~ tol=%e\n',tol);
save('testAXB-v73.mat','A','X','B','-v7.3');
save('testAXB-v7.mat','A','X','B','-v7');
save('testAXB-v6.mat','A','X','B','-v6');
save('testAXB-v4.mat','A','X','B','-v4');
