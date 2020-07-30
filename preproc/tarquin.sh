/Applications/tarquin/tarquingui.app/Contents/MacOS/tarquin \
--input ave_diff_lcm --format lcm \
--lipid_filter false --auto_phase true --auto_ref true \
--crlb_optim false --ref_signals 1h_naa --fs 2000 \
--ft 123260173 --echo .080  --pul_seq mega_press --gnuplot /usr/local/bin/gnuplot \
--output_txt tarquin.txt --output_pdf tarquin.pdf --output_image tarquin_img.pdf \
--ext_pdf true --output_xml tarquin.xml --int_basis megapress_gaba \
--svs_only true --stack_pdf true --te1 .04
