function [ap, ml, dv] = dryadToTemplate(AP, ML, DV, tf)
% [ap, ml, dv] = dryadToTemplate(AP, ML, DV, tf) maps Dryad cell coordinates (AP, ML,
% DV columns of loadDryadCells) into a template frame built by
% makeAnatomyTemplates, using the transform it fitted (T.tf):
%   ap = AP with its spacing scaled by tf.ap_scale about tf.ap_mid (vibratome
%   thickness / nominal 0.1 mm), plus tf.ap_shift; then (ap, DV) rotated tf.theta
%   deg nose down about lambda, then DV + tf.c. ML is unchanged.
%
% Hannah Payne, 2026

ap = tf.ap_mid + (AP - tf.ap_mid) * tf.ap_scale + tf.ap_shift;   % as the sections in buildBirdTemplate
[ap, dv] = rotateNoseDown(ap, DV, tf.theta);
dv = dv + tf.c;
ml = ML;
end
