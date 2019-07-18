#include <iterator>
#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
using namespace sdsl;
using namespace std;
// for (int i = 0; i < m; i++)
// {
//     for (int ci = 0; ci < sigma.length(); ci++)
//     {
//         char c = sigma.at(ci);
//         if (c != query.at(i))
//         {
//             string modified = query.substr(0, i) + c + query.substr(i + 1);
//             occs += count1(fm_index, modified.begin(), modified.end());
//         }
//     }
// }
template <class t_csa, class t_pat_iter>
typename t_csa::size_type count1(
    const t_csa &csa,
    t_pat_iter begin,
    t_pat_iter end,
    csa_tag)
{
    if (end - begin > (typename std::iterator_traits<t_pat_iter>::difference_type)csa.size())
        return 0;

    string sigma = "abcdefghijklmnopqrstuvwxyz";

    string::iterator middle = begin;
    int mod = (end - begin) % 2;
    advance(middle, ((end - begin) / 2) + mod - 1);
    cout << "middle is " << *middle << endl;

    // error is in first half:
    typename t_csa::size_type l_res = 0;
    typename t_csa::size_type r_res = 0;

    typename t_csa::size_type total = 0;

    typename t_csa::size_type result = backward_search(csa, 0, csa.size() - 1, middle, end, l_res, r_res);
    cout << "right side matches " << result << endl;

    if (result)
    {
        for (string::iterator i = middle - 1; i >= begin; i--)
        {
            // the missmatch is at "i"
            typename t_csa::size_type l_res2 = 0;
            typename t_csa::size_type r_res2 = 0;
            int index = (i - begin);
            cout << " char "<< *i << " index = " <<  index << ", middle - 1 = " << (middle - begin - 1) << endl;
            if (i != middle - 1)
            {
                result = backward_search(csa, l_res, r_res, i + 1, middle - 1, l_res2, r_res2);
                cout << "next untill " << *(i + 1) << " results " << result << endl;

                if (!result)
                    continue;
            }
            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                cout << ci << " out of " << sigma.length() << endl;
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res3 = 0;
                typename t_csa::size_type r_res3 = 0;
                char c = sigma.at(ci);
                if (c == *i)
                    continue;                                                     // we do not need to check with the letter in use
                result = backward_search(csa, l_res2, r_res2, c, l_res3, r_res3); // get the range for matches with "c"
                if (!result)
                    continue;
                if (i != begin)
                {                                                                                // only continue if found
                    result = backward_search(csa, l_res3, l_res3, begin, i - 1, l_res3, r_res3); // look for an exact match to the replaced character
                }
                cout << "replace with " << c << " completes " << result << endl;
                total += result;
            }
        }
    }

    // missmatch occours in the right half:

    l_res = 0;
    r_res = 0;

    result = forward_search(csa, 0, csa.size() - 1, begin, middle - 1, l_res, r_res);
    result = 0;
    cout << "left side matches " << result << endl;

    if (result)
    {
        for (string::iterator i = middle; i <= end; i++)
        {
            // the missmatch is at "i"
            typename t_csa::size_type l_res2 = l_res;
            typename t_csa::size_type r_res2 = r_res;
            if (i != middle)
            {
                result = forward_search(csa, l_res, r_res, middle, i - 1, l_res2, r_res2);
                cout << "next untill " << *(i - 1) << " results " << result << endl;

                if (!result)
                    continue;
            }

            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res3 = l_res2;
                typename t_csa::size_type r_res3 = r_res2;
                char c = sigma.at(ci);
                if (c == *i)
                    continue;                                                    // we do not need to check with the letter in use
                result = forward_search(csa, l_res2, r_res2, c, l_res3, r_res3); // get the range for matches with "c"
                if (!result)
                    continue;
                if (i != begin)
                {                                                                             // only continue if found
                    result = forward_search(csa, l_res3, l_res3, i + 1, end, l_res3, r_res3); // look for an exact match from the replaced character
                }
                cout << "replace with " << c << " completes " << result << endl;
                total += result;
            }
        }
    }

    return total;
}

template <class t_csx, class t_pat_iter>
typename t_csx::size_type count1(
    const t_csx &csx,
    t_pat_iter begin,
    t_pat_iter end)
{
    typename t_csx::index_category tag;
    return count1(csx, begin, end, tag);
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage " << argv[0] << " text_file [error_count (0 - 2)] [max_locations] [post_context] [pre_context]" << endl;
        cout << "    This program constructs a very compact FM-index" << endl;
        cout << "    which supports count, locate, and extract queries." << endl;
        cout << "    text_file      Original text file." << endl;
        cout << "    error_count    Exact number of errors in text." << endl;
        cout << "    max_locations  Maximal number of location to report." << endl;
        cout << "    post_context   Maximal length of the reported post-context." << endl;
        cout << "    pre_context    Maximal length of the pre-context." << endl;
        return 1;
    }

    size_t error_count = 0;
    size_t max_locations = 5;
    size_t post_context = 10;
    size_t pre_context = 10;

    if (argc >= 3)
    {
        error_count = atoi(argv[2]);
    }
    if (argc >= 4)
    {
        max_locations = atoi(argv[3]);
    }
    if (argc >= 5)
    {
        post_context = atoi(argv[4]);
    }
    if (argc >= 6)
    {
        pre_context = atoi(argv[5]);
    }
    string index_suffix = ".fm9";
    string index_file = string(argv[1]) + index_suffix;
    csa_wt<wt_huff<rrr_vector<127>>, 512, 1024> fm_index;

    if (!load_from_file(fm_index, index_file))
    {
        ifstream in(argv[1]);
        if (!in)
        {
            cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl;
            return 1;
        }
        cout << "No index " << index_file << " located. Building index now." << endl;
        construct(fm_index, argv[1], 1);     // generate index
        store_to_file(fm_index, index_file); // save it
    }
    cout << "Index construction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB." << endl;
    cout << "Input search terms and press Ctrl-D to exit." << endl;
    string prompt = "\e[0;32m>\e[0m ";
    cout << prompt;
    string query;
    while (getline(cin, query))
    {
        size_t m = query.size();
        // todo: here we count the occurrences
        if (error_count == 0)
        {
            size_t occs = count(fm_index, query.begin(), query.end());
            cout << "# of occurrences: " << occs << endl;
            if (occs > 0)
            {
                cout << "Location and context of first occurrences: " << endl;
                // todo: here we locate
                auto locations = locate(fm_index, query.begin(), query.begin() + m);
                sort(locations.begin(), locations.end());
                for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i)
                {
                    cout << setw(8) << locations[i] << ": ";
                    if (pre_extract > locations[i])
                    {
                        pre_extract = locations[i];
                    }
                    if (locations[i] + m + post_extract > fm_index.size())
                    {
                        post_extract = fm_index.size() - locations[i] - m;
                    }
                    auto s = extract(fm_index, locations[i] - pre_extract, locations[i] + m + post_extract - 1);
                    string pre = s.substr(0, pre_extract);
                    s = s.substr(pre_extract);
                    if (pre.find_last_of('\n') != string::npos)
                    {
                        pre = pre.substr(pre.find_last_of('\n') + 1);
                    }
                    cout << pre;
                    cout << "\e[1;31m";
                    cout << s.substr(0, m);
                    cout << "\e[0m";
                    string context = s.substr(m);
                    cout << context.substr(0, context.find_first_of('\n')) << endl;
                }
            }
        }
        else if (error_count == 1)
        {
            string sigma = "abcdefghijklmnopqrstvuwxyz"; //todo: find a way to achive # of abb
            size_t occs = count1(fm_index, query.begin(), query.end());

            cout << "# of occurrences: " << occs << endl;
        }
        cout << prompt;
    }
    cout << endl;
}
