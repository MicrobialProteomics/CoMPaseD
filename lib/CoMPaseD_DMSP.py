import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
from numpy import mean, hstack, array
from pandas import DataFrame
from os import environ
environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
# import tensorflow as tf
from tensorflow import keras
environ['TF_CPP_MIN_LOG_LEVEL'] = '0'

def run_DeepMSPep(peptides_list, model_path):

    # peptide_list_to_codify = ["PEPTIDE", "PPPEEEPPTIDE", "PEPTIIDE"]
    # max_aa = 81
    def load_pep_and_codify(peptide_list_to_codify, max_aa = 81):
        aa_dict = {'A':1,'R':2,'N':3,'D':4,'C':5,'Q':6,'E':7,'G':8,'H':9,'I':10,'L':11,'K':12,'M':13,'F':14,
                   'P':15,'O':16,'S':17,'U':18,'T':19,'W':20,'Y':21,'V':22}

        pep_codes = list()
        long_pep_counter = 0
        newLines = list()
        # added for error handling of special aa's
        peptides_with_special_aa = list()

        for pep in peptide_list_to_codify:
            if not len(pep) > max_aa:
                current_pep= list()
                for aa in pep:
                    # added error handling
                    try:
                        current_pep.append(aa_dict[aa])
                    except KeyError:
                        peptides_with_special_aa.append(pep)
                pep_codes.append(current_pep)
                newLines.extend([pep])
            else:
                long_pep_counter += 1
        predict_data = keras.preprocessing.sequence.pad_sequences(pep_codes, value=0, padding='post', maxlen=max_aa)
        return predict_data, long_pep_counter, newLines, peptides_with_special_aa

    print(f"---------------------------------------------------------------------------", flush=True)
    print(f"Using DeepMSPeptide to predict peptide detectability", flush=True)
    print(f"\tfor a detailed description please see:", flush=True)
    print(f"\tGuillermo Serrano, Elizabeth Guruceaga and Victor Segura", flush=True)
    print(f"\tDeepMSPeptide: peptide detectability prediction using deep learning.")
    print(f"\tBioinformatics. 36(4):1279-1280, 2019.", flush=True)
    print(f"\tDOI: 10.1093/bioinformatics/btz708", flush=True)
    print("---------------------------------------------------------------------------", flush=True)


    print(f'\tLoading prediction model', flush=True)
    # updated to keras 3; 2024-05-21
    # ToDo: Check if older python versions (3.5-3.11) will work properly with keras 3 / tensorflow 2.16
    model_2_1D = keras.models.load_model(model_path, compile=False)

    print(f'\tLoading input peptides', flush=True)
    predict_data, skipped,  lines, special_peptides = load_pep_and_codify(peptide_list_to_codify = peptides_list, max_aa = model_2_1D.input_shape[1])
    print('\tSuccesfully loaded {0} peptides and skipped {1}'.format(len(lines), str(skipped)), flush=True)

    # in model.predict use verbose=2 to suppress progress bar send to stdout as is appended to the gui progress win
    model_2_1D_pred = model_2_1D.predict(predict_data, verbose=2)
    model_2_1D_pred = hstack((array(lines).reshape(len(lines), 1),model_2_1D_pred)).tolist()

    # export only predictions as True/False option is not required
    Pep_seq = list()
    Pred_output = list()

    # prediction must be converted to numeric
    for pred in model_2_1D_pred:
        Pep_seq.append(pred[0])
        Pred_output.append(float(pred[1]))

    # in case, there were peptide(s) with special characters, set them to the mean prediction
    mean_prediction = mean(Pred_output)

    for special_pep in special_peptides:
        Pep_seq.append(special_pep)
        Pred_output.append(mean_prediction)


    # transform prediction to high values for 'good' peptides
    Pred_output_2 = list()
    for pep_prediction in Pred_output:
        Pred_output_2.append(1.0-float(pep_prediction))

    return_df = DataFrame()
    return_df['peptide'] = Pep_seq
    return_df['DeepMSPep_prediction'] = Pred_output_2

    return return_df
