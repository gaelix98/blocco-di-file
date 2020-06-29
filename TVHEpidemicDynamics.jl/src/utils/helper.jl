"""
    A bunch of helper functions.
"""


"""
    interval_id -> (start_date, end_date)
"""
function evaluateintervals!(mintime, maxtime, intervals, nintervals, Δₘ)
    for i=1:nintervals
        offset = mintime + Δₘ > maxtime ?  maxtime + Dates.Millisecond(1) :  mintime + Δₘ
        #get!(intervals, i, Pair{DateTime, DateTime}(convert(Dates.DateTime, (mintime-delta)), convert(Dates.DateTime, (offset+delta))))
        get!(intervals, i, Pair{DateTime, DateTime}(convert(Dates.DateTime, (mintime)), convert(Dates.DateTime, (offset))))

        mintime = offset
    end
end


"""
 Evaluate the checkin density per time interval,
 counting the number of checkins within an interval.
"""
function evaluatedensity!(intervals, checkinperinterval, df)
    for interval in intervals
        currdf = filter(r -> ((r.UTCtime >= (interval.second.first)) && (r.UTCtime < (interval.second.second))), df)
        push!(
            checkinperinterval,
            interval.first => size(currdf)[1])
    end
end



"""
 Evaluate checkin distribution within a given interval.
"""
function evaluatedistance!(intervals, distancewithinintervals, df)
    for interval in intervals
        currdf = filter(r -> ((r.UTCtime >= (interval.second.first)) && (r.UTCtime < (interval.second.second))), df)
        #println("$(interval) --> $(size(currdf)[1])")

        get!(distancewithinintervals, interval.first, Array{Int, 1}())

        for r₁ in eachrow(currdf)
            for r₂ in eachrow(currdf)
                if getfield(r₁, :row) != getfield(r₂, :row)
                    if r₁.venueid == r₂.venueid
                        push!(
                            get!(distancewithinintervals, interval.first, Array{Int, 1}()),
                            convert(Dates.Second, abs(r₁.UTCtime - r₂.UTCtime)).value #convert(Dates.Second, (mintime))
                            #round(dates[index+1] - dates[index], Dates.Minute).value
                        )
                    end
                end
            end
        end
    end
end


"""
 Evaluate direct contact distribution within a given interval.
"""
function evaluate_direct_contacts_distribution!(intervals, distancewithinintervals, df, δ)
    for interval in intervals
        currdf = filter(r -> ((r.UTCtime >= (interval.second.first)) && (r.UTCtime < (interval.second.second))), df)

        #println("$(interval) --> $(size(currdf)[1])")

        get!(distancewithinintervals, interval.first, Array{Int, 1}())

        users = unique(currdf, :userid)

        for u in eachrow(users)
            uid = u.userid
            ucontacts = filter(r -> r.userid == uid, currdf)

            for uc in eachrow(ucontacts)
                for r in eachrow(currdf)
                    if uc.userid != r.userid
                        if uc.venueid == r.venueid
                            if abs(uc.UTCtime - r.UTCtime) <= δ
                                push!(
                                    get!(distancewithinintervals, interval.first, Array{Int, 1}()),
                                    convert(Dates.Second, abs(uc.UTCtime - r.UTCtime)).value
                                )
                            end
                        end
                    end
                end
            end
        end
    end
end


"""
 Evaluate indirect contact distribution within a given interval.
"""
function evaluate_location_distribution!(intervals, distancewithinintervals, df, δ)
    for interval in intervals
        currdf = filter(r -> ((r.UTCtime >= (interval.second.first)) && (r.UTCtime < (interval.second.second))), df)

        #println("$(interval) --> $(size(currdf)[1])")

        get!(distancewithinintervals, interval.first, Array{Int, 1}())

        users = unique(currdf, :userid)

        for u in eachrow(users)
            uid = u.userid
            ucontacts = filter(r -> r.userid == uid, currdf)

            push!(
                get!(distancewithinintervals, interval.first, Array{Int, 1}()),
                size(unique(ucontacts, :venueid))[1]
            )
        end
    end
end
