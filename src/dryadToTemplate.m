function [ap, ml, dv] = dryadToTemplate(AP, ML, DV, tf)
% [ap, ml, dv] = dryadToTemplate(AP, ML, DV, tf) maps Dryad cells (AP, ML, and DV
% as depth below their section's top, from loadDryadCells with a template's P)
% into the template frame, using the placement makeAnatomyTemplates fitted (T.tf):
%   ap = AP with its spacing scaled by tf.ap_scale about tf.ap_mid (vibratome
%   thickness / nominal 0.1 mm), plus tf.ap_shift, as the template's sections
%   dv = depth + the DV where a section top goes at that ap (tf.top_ap, tf.top_dv)
%   ML is unchanged. Cells beyond tf.top_ap get NaN DV.
%
% Hannah Payne, 2026

ap = tf.ap_mid + (AP - tf.ap_mid) * tf.ap_scale + tf.ap_shift;   % as the sections in buildBirdTemplate
dv = DV + interp1(tf.top_ap, tf.top_dv, ap);
ml = ML;
end
