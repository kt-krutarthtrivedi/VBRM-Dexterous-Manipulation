%   VBRM - Dexterous In-hand Manipulation 
%   Author: Krutarth Ambarish Trivedi | ktrivedi@wpi.edu

clc;
close all;
clear;

%----- Parameters from the given model ------%
time = 2;
rotational_axis = [0;0;0];
theta_dot = 0;
angle_contact_1 = 90;
angle_contact_2 = -90;

%----- Calculating the frames -----%

O_N = [0; 0; 0.05];
T_N = [0; 0; 0.035];

C1_N = [-0.015; 0; 0.05];
C2_N = [0.015; 0; 0.05];

J1_N = [-0.015; 0; 0];
J2_N = [-0.015-(0.05*cos(deg2rad(30))); 0; 0.025];
J3_N = [0.015; 0;0];
J4_N = [0.015+(0.05*cos(deg2rad(30))); 0; 0.025];

%----- Object twist -----%

object_V_N = zeros(3,1);

for i = 1:3
    object_V_N(i,1) = (T_N(i) - O_N(i))/time; 
end

object_W_N = theta_dot .* rotational_axis;

object_twist_N =[object_V_N; object_W_N];

disp("Object twist w.r.t N = ");
disp(object_twist_N);

%----- For contact point 1 -----%

skew_vector_1= C1_N - O_N;
skew_1 = create_skew_matrix(skew_vector_1);
P_1 = [eye(3), transpose(skew_1) ; 
    zeros(3), eye(3)];
disp("P1 = ");
disp(P_1);

contact_1_twist_N = P_1 * object_twist_N;
disp("Contact point 1 twist w.r.t N = ");
disp(contact_1_twist_N);

rot_contact_1_N = find_rotation_y(angle_contact_1);
rot_contact_1_N_bar = [rot_contact_1_N, zeros(3);
    zeros(3), rot_contact_1_N];
disp("Rotation bar of contact 1 w.r.t N = ");
disp(rot_contact_1_N_bar);

rot_N_contact_1_bar = transpose(rot_contact_1_N_bar);
disp("Rotation bar of N w.r.t contact 1 = ");
disp(rot_N_contact_1_bar);

contact_1_twist_contact_1 = rot_N_contact_1_bar * contact_1_twist_N;
disp("Contact point 1 twist w.r.t Contact point 1 = ");
disp(contact_1_twist_contact_1);

%----- For contact point 2 -----%

skew_vector_2= C2_N - O_N;
skew_2 = create_skew_matrix(skew_vector_2);
P_2 = [eye(3), transpose(skew_2) ; 
    zeros(3), eye(3)];
disp("P2 = ");
disp(P_2);

contact_2_twist_N =  P_2 * object_twist_N;
disp("Contact point 2 twist w.r.t N = ");
disp(contact_2_twist_N);

rot_contact_2_N = find_rotation_y(angle_contact_2);
rot_contact_2_N_bar = [rot_contact_2_N, zeros(3);
    zeros(3), rot_contact_2_N];
disp("Rotation bar of contact 2 w.r.t N = ");
disp(rot_contact_2_N_bar);
rot_N_contact_2_bar = transpose(rot_contact_2_N_bar);
disp("Rotation bar of N w.r.t contact 2 = ");
disp(rot_N_contact_2_bar);

contact_2_twist_contact_2 = rot_N_contact_2_bar * contact_2_twist_N;
disp("Contact point 2 twist w.r.t Contact point 2 = ");
disp(contact_2_twist_contact_2);

%----- Grasp Matrix -----%

Grasp_Matrix = [rot_N_contact_1_bar * P_1;
    rot_N_contact_2_bar * P_2];
disp("Grasp Matrix = ");
disp(Grasp_Matrix);

contacts_twist_contacts = Grasp_Matrix * object_twist_N;
disp("Contact twists w.r.t their respective contact frames = ")
disp(contacts_twist_contacts);

%---------- Hand Jacobian -------------%

C1_N = [-0.015; 0; 0.05];
C2_N = [0.015; 0; 0.05];

J1_N = [-0.015; 0; 0];
J2_N = [-0.015-(0.05*cos(deg2rad(30))); 0; 0.025];
J3_N = [0.015; 0;0];
J4_N = [0.015+(0.05*cos(deg2rad(30))); 0; 0.025];

W_joints = [0;1;0];
W_J1 = W_joints;
W_J2 = W_joints;
W_J3 = W_joints;
W_J4 = W_joints;

skew_W_J1 = create_skew_matrix(W_J1);
skew_W_J2 = create_skew_matrix(W_J2);
skew_W_J3 = create_skew_matrix(W_J3);
skew_W_J4 = create_skew_matrix(W_J4);

V_J1 = skew_W_J1 * (C1_N - J1_N);
V_J2 = skew_W_J2 * (C1_N - J2_N);
V_J3 = skew_W_J3 * (C2_N - J3_N);
V_J4 = skew_W_J4 * (C2_N - J4_N);

Jacobian_contact_1 = [V_J1, V_J2;
    W_J1, W_J2];

Jacobian_contact_2 = [V_J3, V_J4;
    W_J3, W_J4];

disp("Jacobian Matrix for Contact Point 1 = ");
disp(Jacobian_contact_1);
disp("Jacobian Matrix for Contact Point 2 = ");
disp(Jacobian_contact_2);

Hand_Jacobian = blkdiag(rot_N_contact_1_bar * Jacobian_contact_1, rot_N_contact_2_bar * Jacobian_contact_2);

disp("Hand Jacobian = ");
disp(Hand_Jacobian);

%----------------- Jacobian Hand_Object --------------%

Hand_Object_Jacobian = pinv(Grasp_Matrix) * Hand_Jacobian;

disp("Hand Object Jacobian = ");
disp(Hand_Object_Jacobian);

%------------- Finding Joint angles ----------%

joint_angles = pinv(Hand_Jacobian) * Grasp_Matrix * object_twist_N;
disp("Joint Angles = ");
disp(joint_angles);

%-------------- User defined functions --------------%

% calculate a skew symmetric matrix from a 3d vector
function S = create_skew_matrix(vector)
    S(1,2) = -vector(3);
    S(1,3) = vector(2);
    S(2,3) = -vector(1);
    S(2,1) = vector(3);
    S(3,1) = -vector(2);
    S(3,2) = vector(1);
end

% calculate the rotation around Y axis for the given angle
function R = find_rotation_y(theta)
    R = [cos(deg2rad(theta)), 0, sin(deg2rad(theta));
        0, 1, 0;
        -sin(deg2rad(theta)),0, cos(deg2rad(theta))];
end