import main as preproc
import compile_biokg
import build_benchmarks

def main():
    print('Preprocessing data')
    preproc.preprocess_graph()
    print('Compiling graph')
    compile_biokg.compile_graph()
    print('Building benchmarks')
    build_benchmarks.main()


if __name__ == '__main__':
    main()