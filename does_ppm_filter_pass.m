function exp = does_ppm_filter_pass(exp)

max_abs_ppm_err = max(abs(exp.data.ppm_err_TMTc.*exp.binary_abundand_Ynmin1_positions(:,1:12)),[],2);
exp.passed_ppm_filter = max_abs_ppm_err < 10;

end