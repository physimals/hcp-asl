from nbclient import execute
import nbformat
from nbparameterise import extract_parameters, parameter_values, replace_definitions

from pathlib import Path

def create_qc_report(subject_dir, outdir):
    # get location of this file
    file_name = Path(__file__).resolve(strict=True)

    # derive location of template using location of this file
    hcpasl_dir = file_name.parent
    template_nb = nbformat.read(str(hcpasl_dir/"report_template.ipynb"), as_version=4)

    # extract original parameters from the template notebook
    orig_parameters = extract_parameters(template_nb)

    # replace with parameters for this subject and execute notebook
    new_parameters = parameter_values(orig_parameters,
                                      subject_dir=str(subject_dir),
                                      outdir=str(outdir))
    new_nb = replace_definitions(template_nb, new_parameters, execute=False)
    _ = execute(new_nb)

    # save notebook in subject's main output directory
    new_nb_name = Path(subject_dir)/outdir/"report.ipynb"
    with open(new_nb_name, "w") as f:
        nbformat.write(new_nb, f)
        