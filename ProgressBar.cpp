// The MIT License (MIT)
//
// Copyright (c) 2019 Luigi Pertoldi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// ============================================================================
//  ___   ___   ___   __    ___   ____  __   __   ___    __    ___
// | |_) | |_) / / \ / /`_ | |_) | |_  ( (` ( (` | |_)  / /\  | |_)
// |_|   |_| \ \_\_/ \_\_/ |_| \ |_|__ _)_) _)_) |_|_) /_/--\ |_| \_
//
// Very simple progress bar for c++ loops with internal running variable
//
// Author: Luigi Pertoldi
// Created: 3 dic 2016
//
// Notes: The bar must be used when there's no other possible source of output
//        inside the for loop
//

#include "ProgressBar.h"

ProgressBar::ProgressBar() :
    progress(0),
    n_cycles(0),
    last_perc(0),
    do_show_bar(true),
    update_is_called(false),
    done_char("#"),
    todo_char(" "),
    opening_bracket_char("["),
    closing_bracket_char("]") {}

ProgressBar::ProgressBar(int n, bool showbar) :
    progress(0),
    n_cycles(n),
    last_perc(0),
    do_show_bar(showbar),
    update_is_called(false),
    done_char("#"),
    todo_char(" "),
    opening_bracket_char("["),
    closing_bracket_char("]") {}

void ProgressBar::reset() {
    progress = 0,
    update_is_called = false;
    last_perc = 0;
    return;
}

void ProgressBar::set_niter(int niter) {
    if (niter <= 0) throw std::invalid_argument(
        "ProgressBar::set_niter: number of iterations null or negative");
    n_cycles = niter;
    return;
}

void ProgressBar::update() {

    if (n_cycles == 0) throw std::runtime_error(
            "ProgressBar::update: number of cycles not set");

    if (!update_is_called) {
        if (do_show_bar == true) {
            std::cout << opening_bracket_char;
            for (int _ = 0; _ < 50; _++) std::cout << todo_char;
            std::cout << closing_bracket_char << " 0%";
        }
        else std::cout << "0%";
    }
    update_is_called = true;

    int perc = 0;

    // compute percentage, if did not change, do nothing and return
    perc = progress*100./(n_cycles-1);
    if (perc < last_perc) return;

    // update percentage each unit
    if (perc == last_perc + 1) {
        // erase the correct  number of characters
        if      (perc <= 10)                std::cout << "\b\b"   << perc << '%';
        else if (perc  > 10 and perc < 100) std::cout << "\b\b\b" << perc << '%';
        else if (perc == 100)               std::cout << "\b\b\b" << perc << '%';
    }
    if (do_show_bar == true) {
        // update bar every ten units
        if (perc % 2 == 0) {
            // erase closing bracket
            std::cout << std::string(closing_bracket_char.size(), '\b');
            // erase trailing percentage characters
            if      (perc  < 10)               std::cout << "\b\b\b";
            else if (perc >= 10 && perc < 100) std::cout << "\b\b\b\b";
            else if (perc == 100)              std::cout << "\b\b\b\b\b";

            // erase 'todo_char'
            for (int j = 0; j < 50-(perc-1)/2; ++j) {
                std::cout << std::string(todo_char.size(), '\b');
            }

            // add one additional 'done_char'
            if (perc == 0) std::cout << todo_char;
            else           std::cout << done_char;

            // refill with 'todo_char'
            for (int j = 0; j < 50-(perc-1)/2-1; ++j) std::cout << todo_char;

            // readd trailing percentage characters
            std::cout << closing_bracket_char << ' ' << perc << '%';
        }
    }
    last_perc = perc;
    ++progress;
    std::cout << std::flush;

    return;
}
