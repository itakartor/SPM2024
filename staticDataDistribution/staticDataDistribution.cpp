
// cyclic
auto cyclic = [&] (const index_t& id) -> void {
    for (index_t row = id; row < m; row += num_threads) {
        value_t accum = value_t(0);
        for (index_t col = 0; col < n; col++) {
            accum += A[row*n+col]*x[col];
        }
        b[row] = accum;
    }
};